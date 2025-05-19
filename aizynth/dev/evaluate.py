# python

import json
import logging

from rdkit import Chem
from rdkit.Chem.rdChemReactions import ReactionFromSmarts

from aizynth.classifier.classifier import transform_retro_to_forward

# Create logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# Add console handler
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


def add_to_label_dict(uspto_templates: list[dict], threshold: float, label_dict: dict, rxn_name: str,
                      rxn_smirks: str, parent_name: str | None = None, parent_smirks: str | None = None) -> dict:
    """Add the uspto templates to the label dict for the given threshold."""
    for uspto_temp in uspto_templates:
        temp_code: int = uspto_temp['template_code']
        rxn_sim: float = uspto_temp['rxn_sim']
        if rxn_sim >= threshold:
            to_add: dict = {
                'sim': rxn_sim,
                'name': rxn_name,
                'smirks': rxn_smirks
            }
            if parent_name is not None:
                to_add['parent_name'] = parent_name
                to_add['parent_smirks'] = parent_smirks
            if temp_code not in label_dict:
                label_dict[temp_code] = [to_add]
            else:
                label_dict[temp_code].append(to_add)
    return label_dict


def get_label_dict(labels: list[dict], test: list[dict], threshold: float, method: str) -> dict:
    """Produce a dictionary of labels for the given threshold."""
    label_dict = {}
    if method == 'all':
        logger.info(f"Formatting all labels for {threshold}")
        for group in labels:
            # parent
            parent_name = group['name']
            parent_smirks = group['smirks']
            label_dict = add_to_label_dict(group['uspto_templates'], threshold, label_dict, parent_name, parent_smirks)
            # children
            for child in group['matched']:
                child_name = child['name']
                child_smirks = child['smirks']
                label_dict = add_to_label_dict(child['uspto_templates'], threshold, label_dict, child_name,
                                               child_smirks)
    elif method == 'parent':
        logger.info(f"Formatting using parent and child labels for {threshold}")
        for group in labels:
            # parent
            parent_name = group['name']
            parent_smirks = group['smirks']
            label_dict = add_to_label_dict(group['uspto_templates'], threshold, label_dict, parent_name, parent_smirks,
                                           parent_name, parent_smirks)
            # children
            for child in group['matched']:
                child_name = child['name']
                child_smirks = child['smirks']
                label_dict = add_to_label_dict(child['uspto_templates'], threshold, label_dict, child_name,
                                               child_smirks, parent_name, parent_smirks)
    elif method == 'only_parent':
        logger.info(f"Formatting using only parent labels for {threshold}")
        for group in labels:
            # parent
            parent_name = group['name']
            parent_smirks = group['smirks']
            label_dict = add_to_label_dict(group['uspto_templates'], threshold, label_dict, parent_name, parent_smirks)

    elif method == 'exact':
        logger.info(f"Formatting using exact template matches from test data for {threshold}")
        for rxn in test:
            if 'route_metadata' in rxn and 'template_code' in rxn['route_metadata']:
                temp_code = rxn['route_metadata']['template_code']
                smirks = rxn['route_metadata'].get('template', None)
                if smirks:
                    if temp_code not in label_dict:
                        label_dict[temp_code] = [{
                            'sim': 1.0,
                            'name': None,
                            'smirks': transform_retro_to_forward(smirks)
                        }]

    # order the labels by decreasing similarity
    logger.info(f"Sorting labels for {threshold}")
    for key in label_dict:
        label_dict[key] = sorted(label_dict[key], key=lambda x: x['sim'], reverse=True)
    return label_dict


def can_be_sanitized(mol: Chem.Mol) -> bool:
    if type(mol) != Chem.Mol:
        logger.error(f"Expected a Chem.Mol object, got {type(mol)}.")  # Make sure it's a Chem.Mol object
    try:
        Chem.SanitizeMol(mol)
        return True
    except:
        return False


def evaluate_product_list(products: list[tuple[Chem.Mol]], product_inchi: str, label: dict) -> bool:
    """Evaluate the product list against the product."""
    found_product = False
    if len(products) == 0:
        logger.info(f"No products found for {label['name']}")
        return found_product
    if len(products) > 1 or len(products[0]) > 1:
        logger.warning(f"More than one product found for {label['name']}")
        pred_prods = [product[0] for product in products if can_be_sanitized(product[0])]
        if len(pred_prods) == 0:
            logger.info(f"No sanitizeable products found for {label['name']}")
            return found_product
    elif len(products) == 1 and len(products[0]) == 1:
        logger.info(f"One product found for {label['name']}")
        pred_prods = [products[0][0]]
    else:
        logger.warning(f"Unexpected number of products found for {label['name']}")
        return found_product

    for product in pred_prods:
        pred_prod_inchi: str = Chem.MolToInchiKey(product)
        if pred_prod_inchi == product_inchi:
            found_product = True
            break
    return found_product


def pred_product(product: str, reactants: list[str], sim_labels: list[dict], method: str, true_label: str) -> dict:
    """Try to predict the product from the reactants using templates."""
    found_product = False
    rxn_name_worked: tuple[str] | str | None = None
    n_tried = 0
    rxn_sim_worked: float | None = None
    smirks: tuple[str] | str | None = None
    true_label_accessible: bool = len(true_label) > 0

    product_inchi: str = Chem.MolToInchiKey(Chem.MolFromSmiles(product))

    # if len(reactants) != 2:
    #     logger.error(f"Expected two reactants, got {len(reactants)}.")
    #     results = {
    #         'found_product': found_product,
    #         'n_tried': n_tried,
    #         'rxn_name_worked': None,
    #         'rxn_sim_worked': None,
    #         'smirks_worked': None,
    #         'true_label_accessible': true_label_accessible
    #     }
    #     return results

    for label in sim_labels:
        n_tried += 1
        if label['name'] == 'Suzuki_coupling':
            pass
        reaction: Chem.rdChemReactions = ReactionFromSmarts(label['smirks'])
        if len(reaction.GetReactants()) != len(reactants):
            logger.error(f"Expected {len(reactants)} reactants, got {len(reaction.GetReactants())}.")
            continue
        if len(reactants) == 1:
            r1 = Chem.MolFromSmiles(reactants[0])
            r2 = None
            products = reaction.RunReactants((r1,))
            found_product1 = evaluate_product_list(products, product_inchi, label)
            found_product2 = False
        elif len(reactants) == 2:
            r1 = Chem.MolFromSmiles(reactants[0])
            r2 = Chem.MolFromSmiles(reactants[1])
            try:
                products1: list[Chem.Mol] = reaction.RunReactants((r1, r2))
                products2: list[Chem.Mol] = reaction.RunReactants((r2, r1))
            except ValueError:
                logger.error(f"Error running reaction {label['name']} with reactants {reactants}.")
                continue

            found_product1 = evaluate_product_list(products1, product_inchi, label)
            found_product2 = evaluate_product_list(products2, product_inchi, label)
        else:
            logger.error(f"Unexpected number of reactants for {label['name']}")
            continue
        if method == 'parent':
            if found_product1 or found_product2:
                # additionally evaluate the products with the parent label
                parent_reaction: Chem.rdChemReactions = ReactionFromSmarts(label['parent_smirks'])
                if len(parent_reaction.GetReactants()) != len(reactants):
                    logger.error(
                        f"Expected {len(reactants)} parent reactants, got {len(parent_reaction.GetReactants())}.")
                    continue
                if len(reactants) == 1:
                    parent_products1 = parent_reaction.RunReactants((r1,))
                    found_product1 = evaluate_product_list(parent_products1, product_inchi, label)
                    found_product2 = False
                elif len(reactants) == 2:
                    parent_products1: list[Chem.Mol] = parent_reaction.RunReactants((r1, r2))
                    parent_products2: list[Chem.Mol] = parent_reaction.RunReactants((r2, r1))
                    found_product1 = evaluate_product_list(parent_products1, product_inchi, label)
                    found_product2 = evaluate_product_list(parent_products2, product_inchi, label)
                else:
                    logger.error(f"Unexpected number of parent reactants for {label['name']}")
                    continue
                if found_product1 or found_product2:
                    found_product = True
                    rxn_name_worked = (label['name'], label['parent_name'])
                    rxn_sim_worked = label['sim']
                    smirks = (label['smirks'], label['parent_smirks'])
                    break
        else:
            if found_product1 or found_product2:
                found_product = True
                rxn_name_worked = label['name']
                rxn_sim_worked = label['sim']
                smirks = label['smirks']
                break
    results = {
        'found_product': found_product,
        'n_tried': n_tried,
        'rxn_name_worked': rxn_name_worked,
        'rxn_sim_worked': rxn_sim_worked,
        'smirks_worked': smirks,
        'true_label_accessible': true_label_accessible
    }
    return results


def evaluate_rxn(rxn: dict, sim_labels: list[dict], method: str) -> dict:
    """Evaluate the reaction against the labels."""
    product: str = rxn['product']
    reactants: list[str] = rxn['reactants']
    true_label: str = rxn['label']
    return pred_product(product, reactants, sim_labels, method, true_label)


def evaluate(test: list[dict], orig_labels: list[dict]) -> list[dict]:
    """Evaluate the classifier on the test data."""
    overall_results = []
    for threshold in [0.2, 0.25, 0.3, 0.3125, 0.35, 0.4]:  # get the labels for the given threshold
        for method in ['exact', 'all', 'parent', 'only_parent']:
            logger.info(f"Evaluating with threshold {threshold} and method {method}")
            per_rxn_results = []
            labels: dict = get_label_dict(orig_labels, test, threshold, method)
            for i, rxn in enumerate(test):
                if 'route_metadata' not in rxn:
                    logger.warning(f"Skipping reaction {rxn} because it does not have route metadata.")
                    continue
                template_code: int = rxn['route_metadata']['template_code']
                if template_code not in labels:
                    logger.warning(f"Template code {template_code} does not have any similar smirks.")
                    results = {
                        'found_product': False,
                        'n_tried': 0,
                        'rxn_name_worked': None,
                        'rxn_sim_worked': None,
                        'smirks_worked': None,
                        'true_label_accessible': None
                    }
                    per_rxn_results.append(results)
                    continue
                sim_labels: list[dict] = labels[template_code]
                if i == 33 and method == 'all':
                    pass
                results: dict = evaluate_rxn(rxn, sim_labels, method)
                results['rxn_i'] = i
                per_rxn_results.append(results)
                setup_results = {
                    'threshold': threshold,
                    'method': method,
                    'n_rxns_evaluated': len(per_rxn_results),
                    'n_products_found_correctly': sum(
                        [result['found_product'] for result in per_rxn_results if result['found_product'] is not None]),
                    'n_products_not_found': sum([not result['found_product'] for result in per_rxn_results if
                                                 result['found_product'] is not None]),
                    'n_products_found_but_label_inaccessible': sum(
                        [result['found_product'] and not result['true_label_accessible'] for result in per_rxn_results
                         if result['found_product'] is not None and result['true_label_accessible'] is not None]),
                    'n_products_not_found_but_label_accessible': sum(
                        [not result['found_product'] and result['true_label_accessible'] for result in per_rxn_results
                         if result['found_product'] is not None and result['true_label_accessible'] is not None]),
                    'n_products_not_found_and_label_inaccessible': sum(
                        [not result['found_product'] and not result['true_label_accessible'] for result in
                         per_rxn_results if
                         result['found_product'] is not None and result['true_label_accessible'] is not None]),
                    'n_products_found_and_label_accessible': sum(
                        [result['found_product'] and result['true_label_accessible'] for result in per_rxn_results if
                         result['found_product'] is not None and result['true_label_accessible'] is not None]),
                    'n_overall_label_accessible': sum(
                        [result['true_label_accessible'] for result in per_rxn_results if
                         result['true_label_accessible'] is not None]),
                    'avg_attempts_per_reaction': sum(
                        [result['n_tried'] for result in per_rxn_results if
                         result['n_tried'] is not None and result['n_tried'] > 0]) / len(per_rxn_results),
                    'per_reaction_results': per_rxn_results
                }
            overall_results.append(setup_results)
    return overall_results


def main():
    test_path = '../data/2025-04-04_rxn_classifer_91_test_data_labeled.json'
    for label_path in ['../data/dev/grouped_smirks_with_codes_concat_True_fp_morgan_threshold_0.2.json',
                       '../data/dev/grouped_smirks_with_codes_concat_True_fp_maccs_threshold_0.2.json']:
        with open(label_path, 'r') as f:
            labels = json.load(f)

        with open(test_path, 'r') as f:
            test = json.load(f)

        results = evaluate(test, labels)

        with open(label_path.replace('.json', '_results.json'), 'w') as f:
            json.dump(results, f, indent=4)


if __name__ == '__main__':
    main()
