# python

import itertools
import json
import logging
import os
from collections import Counter
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from scipy.spatial.distance import jaccard
from sklearn.metrics.pairwise import cosine_similarity


def setup_logging(log_level=logging.INFO, log_file=None) -> logging.Logger:
    """Set up logging configuration."""
    # Create logger
    logger = logging.getLogger()
    logger.setLevel(log_level)

    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # Add console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # Add file handler if log_file is specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    return logger


logger = setup_logging()


def calc_cosine_similarity(a: np.array, b: np.array) -> float:
    """Calculate cosine similarity between two vectors"""
    # Reshape for sklearn format
    a_reshaped = a.reshape(1, -1)
    b_reshaped = b.reshape(1, -1)
    return cosine_similarity(a_reshaped, b_reshaped)[0][0]


def calc_jaccard_similarity(a: np.array, b: np.array) -> float:
    """Calculate taniomoto similarity between two vectors"""
    return 1 - jaccard(a, b)


def morgan_fp(mol: Chem.rdchem.Mol) -> npt.NDArray[Any]:
    """Get morgan fingerprint"""
    fp1 = AllChem.GetMorganFingerprintAsBitVect(
        mol, useChirality=True, radius=2, nBits=1024
    )
    vec1 = np.array(fp1)
    return vec1


def maccs_fp(mol: Chem.rdchem.Mol) -> npt.NDArray[Any]:
    """Get MACCS fingerprint"""
    return np.array(MACCSkeys.GenMACCSKeys(mol))


def get_fp(reaction: str, fp: str = "MACCS", concatenate: bool = True) -> npt.NDArray[Any]:
    reactant_str, product_str = reaction.split(">>")
    reactants = reactant_str.split(".")
    products = product_str.split(".")
    logger.debug(f"Reactants: {reactants}")
    logger.debug(f"Products: {products}")
    reactant_mols = [Chem.MolFromSmarts(reactant) for reactant in reactants]
    product_mols = [Chem.MolFromSmarts(product) for product in products]

    for mol in reactant_mols:
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSSSR(mol)

    for mol in product_mols:
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSSSR(mol)

    # logger.info(f"Reactants: {reactant_mols}")
    # logger.info(f"Products: {product_mols}")

    if fp.lower() == "maccs":
        reactant_fp = np.sum(np.array([maccs_fp(mol) for mol in reactant_mols]), axis=0)
        product_fp = np.sum(np.array([maccs_fp(mol) for mol in product_mols]), axis=0)
    elif fp.lower() == "morgan":
        reactant_fp = np.sum(
            np.array([morgan_fp(mol) for mol in reactant_mols]), axis=0
        )
        product_fp = np.sum(np.array([morgan_fp(mol) for mol in product_mols]), axis=0)
    else:
        raise KeyError(
            f"Fingerprint {fp} is not yet supported. Choose between MACCS and Morgan"
        )

    if concatenate:
        reaction_fp = np.concatenate((reactant_fp, product_fp))
    else:
        reaction_fp = np.sum((reactant_fp, product_fp), axis=0)

    return reaction_fp


def transform_retro_to_forward(reaction: str) -> str:
    """Given a retrosynthesis reaction, transform it into a forward one."""
    parts = reaction.split('>>')
    product, reactants = parts
    return f"{reactants}>>{product}"


def load_data(path: str) -> pd.DataFrame:
    """Load data from .json into pd.DataFrame"""
    assert os.path.exists(path), logger.error(f'{path} does not exist')
    df = pd.read_json(path)
    df.dropna(inplace=True, subset=['reaction'])
    df['reaction'] = df['reaction'].apply(lambda x: transform_retro_to_forward(x))
    logger.info(f'Loaded {len(df)} reactions from {path}')
    return df


def classify_reaction(reaction_fp: np.array, smirks: pd.DataFrame, similarity_metric: str, top_n: int,
                      threshold: float) -> list[dict]:
    """Classify a reaction based on fingerprint similarity. Returns a list of reaction labels."""
    similarities = []

    for _, smirk in smirks.iterrows():
        if similarity_metric == 'cosine':
            sim = calc_cosine_similarity(reaction_fp, smirk['reaction_fp'])
        elif similarity_metric == 'jaccard':
            sim = calc_jaccard_similarity(reaction_fp, smirk['reaction_fp'])
        else:
            raise ValueError(f"Unknown similarity metric: {similarity_metric}")

        similarities.append({
            'name': smirk['name'],
            'similarity': sim
        })

    # Sort by similarity in descending order
    similarities.sort(key=lambda x: x['similarity'], reverse=True)

    # Return top N matches above threshold
    return [s for s in similarities[:top_n] if s['similarity'] >= threshold]


def analyze_by_reaction_type(results: list[dict], fp_type: str, concatenate: bool, similarity_metric: str, top_n: int,
                             threshold: float, overall_accuracy: float,
                             overall_no_label_accuracy: float) -> pd.DataFrame:
    """Analyze performance by reaction type"""
    df = pd.DataFrame(results)
    df['fp_type'] = fp_type
    df['concatenate'] = concatenate
    df['similarity_metric'] = similarity_metric
    df['top_n'] = top_n
    df['threshold'] = threshold

    counts = Counter(df['true_label'])
    accuracy_by_type = df.groupby('true_label')['correct'].mean()
    pred_counts = Counter(list(itertools.chain.from_iterable(df['predicted_labels'])))

    # Plot results
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(25, 8))

    # Frequency plot
    counts_df = pd.DataFrame(list(counts.items()), columns=['reaction_type', 'count'])
    counts_df = counts_df.sort_values('count', ascending=False)
    counts_df['reaction_type'] = counts_df['reaction_type'].astype(str)
    ax1.bar(x=counts_df['reaction_type'], height=counts_df['count'])
    ax1.set_title('Reaction Type Frequency')
    ax1.set_xticklabels(counts_df['reaction_type'], rotation=45, ha='right')

    # Accuracy plot
    accuracy_df = pd.DataFrame(accuracy_by_type).reset_index()
    accuracy_df['true_label'] = accuracy_df['true_label'].astype(str)
    ax2.bar(accuracy_df['true_label'], accuracy_df['correct'])
    ax2.set_title('Accuracy by Reaction Type')
    ax2.set_xticklabels(accuracy_df['true_label'], rotation=45, ha='right')

    # overall frequency of pred labels
    pred_counts_df = pd.DataFrame(list(pred_counts.items()), columns=['predicted_labels', 'count'])
    pred_counts_df = pred_counts_df.sort_values('count', ascending=False)
    pred_counts_df['predicted_labels'] = pred_counts_df['predicted_labels'].astype(str)
    ax3.bar(x=pred_counts_df['predicted_labels'], height=pred_counts_df['count'])
    ax3.set_title('Predicted Labels Frequency')
    ax3.set_xticklabels(pred_counts_df['predicted_labels'], rotation=45, ha='right')

    # overall title
    fig.suptitle(
        f"FP: {fp_type}, Concatenate: {concatenate}, Metric: {similarity_metric}, Top N: {top_n}, Threshold: {threshold}\nAccuracy: {overall_accuracy:.2f}, No Label Accuracy: {overall_no_label_accuracy:.2f}",
        fontsize=16)

    plt.tight_layout()
    plt.savefig(
        f'../data/rxn_type/reaction_type_{fp_type}_{similarity_metric}_{top_n}_{threshold}_concat_{concatenate}.png')
    return df


def evaluate_classifier(test_data: pd.DataFrame, smirks: pd.DataFrame, similarity_metric: str, top_n: int,
                        threshold: float, fp_type: str, concatenate: bool) -> dict:
    """Evaluate classifier w set hyperparameters"""
    results = []
    correct_count = 0
    total_labeled = 0
    no_label_correct_count = 0
    no_label_total = 0

    for _, row in test_data.iterrows():
        total_labeled += 1

        # Get predictions
        pred_results: list[dict] = classify_reaction(
            row['reaction_fp'],
            smirks,
            similarity_metric,
            top_n,
            threshold
        )

        correct = False
        if row['label'] == "":
            no_label_total += 1
            true_label = "no_label"
            if len(pred_results) == 0:
                correct_count += 1  # correctly predicted no reaction
                no_label_correct_count += 1
                pred_labels = ["no_label"]
                correct = True
            else:
                pred_labels = [p['name'] for p in pred_results]
                correct = False
        else:
            # Check if true label is in predicted labels
            pred_labels = [p['name'] for p in pred_results]
            true_label = row['label']

            if true_label in pred_labels:
                correct_count += 1
                correct = True

        # Store detailed result
        results.append({
            'true_label': true_label,
            'predicted_labels': pred_labels,
            'correct': correct,
            'top_similarity': pred_results[0]['similarity'] if pred_results else 0.0
        })

    # Calculate metrics
    accuracy = correct_count / total_labeled if total_labeled > 0 else 0.0
    no_label_accuracy = no_label_correct_count / no_label_total

    if accuracy > 0.4:
        results_by_rxn_type: pd.DataFrame = analyze_by_reaction_type(results, fp_type=fp_type, concatenate=concatenate,
                                                                     similarity_metric=similarity_metric, top_n=top_n,
                                                                     threshold=threshold, overall_accuracy=accuracy,
                                                                     overall_no_label_accuracy=no_label_accuracy)

    return {
        'accuracy': accuracy,
        'correct_count': correct_count,
        'total_labeled': total_labeled,
        'no_label_accuracy': no_label_accuracy,
        'no_label_correct_count': no_label_correct_count,
        'no_label_total': no_label_total,
    }


def evaluate_w_hyperparameters(test_data: pd.DataFrame, smirks: pd.DataFrame, fp_type: str,
                               concatenate: bool) -> tuple[pd.DataFrame, dict]:
    """Evaluate classifier exploring w different hyperparameters"""
    similarity_metrics = ['cosine', 'jaccard']
    top_n_values = [1, 5, 10]
    thresholds = [0.0, 0.3, 0.4, 0.5]
    results = []
    for similarity_metric in similarity_metrics:
        for top_n in top_n_values:
            for threshold in thresholds:
                logger.info(f"Evaluating: {similarity_metric}, top_n={top_n}, threshold={threshold}")
                # Evaluate classifier with current parameters
                eval_result: dict = evaluate_classifier(
                    test_data,
                    smirks,
                    similarity_metric,
                    top_n,
                    threshold,
                    fp_type=fp_type,
                    concatenate=concatenate
                )
                # Store results
                results.append({
                    'similarity_metric': similarity_metric,
                    'top_n': top_n,
                    'threshold': threshold,
                    'accuracy': eval_result['accuracy'],
                    'correct_count': eval_result['correct_count'],
                    'total_labeled': eval_result['total_labeled'],
                    'no_label_accuracy': eval_result['no_label_accuracy'],
                    'no_label_correct_count': eval_result['no_label_correct_count'],
                    'no_label_total': eval_result['no_label_total'],
                })
    results_df = pd.DataFrame(results)
    # Find best configuration
    best_idx = results_df['accuracy'].idxmax()
    best_config = results_df.iloc[best_idx]
    logger.info(f"Best config: {best_config}")
    return results_df, best_config


def plot_overall_results(all_results_df: pd.DataFrame,
                         path_to_save: str = '../data/fingerprint_configurations.png') -> str:
    """Plot all results."""
    plt.figure(figsize=(10, 6))
    bars = plt.bar(
        range(len(all_results_df)),
        all_results_df['accuracy'],
        tick_label=[f"{r['fp_type']}, concat={r['concatenate']}" for _, r in all_results_df.iterrows()]
    )

    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2., height + 0.01,
                 f'{height:.3f}', ha='center', va='bottom')

    plt.xlabel('Fingerprint Configuration')
    plt.ylabel('Best Accuracy')
    plt.title('Best Accuracy by Fingerprint Configuration')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(path_to_save)

    return path_to_save


def get_similarity_heatmap(smirks: pd.DataFrame) -> pd.DataFrame:
    """Calculate similarity heatmap of all smirks"""
    results = []
    sim_metrics = ['cosine', 'jaccard']
    for _, row in smirks.iterrows():
        for _, row_next in smirks.iterrows():
            for similarity_metric in sim_metrics:
                fp1 = row['reaction_fp']
                fp2 = row_next['reaction_fp']
                if similarity_metric == 'cosine':
                    sim = calc_cosine_similarity(fp1, fp2)
                elif similarity_metric == 'jaccard':
                    sim = calc_jaccard_similarity(fp1, fp2)
                result = {
                    'rxn1': row['name'],
                    'rxn2': row_next['name'],
                    'similarity': sim,
                    'similarity_metric': similarity_metric,
                }
                results.append(result)
    return pd.DataFrame(results)


def plot_heatmaps(smirks_similarities: list[pd.DataFrame]) -> list[str]:
    """Plot heatmaps of all smirks"""
    # Ensure output directory exists
    output_dir = '../data/heatmaps'
    os.makedirs(output_dir, exist_ok=True)

    saved_paths = []

    # Process each DataFrame of similarity metrics
    for df in smirks_similarities:
        # Extract configuration info
        fp_type = df['fp_type'].iloc[0]
        concatenate = df['concatenate'].iloc[0]

        # Process each similarity metric separately
        for metric in df['similarity_metric'].unique():
            # Filter data for this metric
            metric_data = df[df['similarity_metric'] == metric]

            # Create a pivot table for the heatmap
            pivot_data = metric_data.pivot(index='rxn1', columns='rxn2', values='similarity')

            # Create the figure
            plt.figure(figsize=(12, 10))

            # Create the heatmap
            sns.heatmap(
                pivot_data,
                annot=False,
                cmap='viridis',  # Color map
                vmin=0,  # Minimum similarity value
                vmax=1,  # Maximum similarity value
                linewidths=0.5,  # Add lines between cells
                fmt=".2f"  # Format for annotations (2 decimal places)
            )

            # Set title and labels
            concat_str = "concatenated" if concatenate else "summed"
            plt.title(
                f'Reaction Similarity Heatmap\n{metric.capitalize()} Similarity - {fp_type} Fingerprints ({concat_str})')
            plt.tight_layout()

            # Save the figure
            filename = f'heatmap_{metric}_{fp_type}_{"concat" if concatenate else "sum"}.png'
            filepath = os.path.join(output_dir, filename)
            plt.savefig(filepath)
            plt.close()

            saved_paths.append(filepath)

    return saved_paths


def run_classifier(test_data: pd.DataFrame, smirks: pd.DataFrame) -> pd.DataFrame:
    """Run evaluation of classifer with different fingerprint configurations"""
    # Configurations to test
    fp_types = ['MACCS', 'morgan']
    concatenate_options = [True, False]
    all_results: list[pd.DataFrame] = []
    smirks_sim_heatmaps: list[pd.DataFrame] = []
    for fp_type in fp_types:
        for concatenate in concatenate_options:
            logger.info(f"Testing configuration: {fp_type}, concatenate={concatenate}")
            test_data['reaction_fp'] = test_data['reaction'].apply(
                lambda r: get_fp(r, fp=fp_type, concatenate=concatenate)
            )
            smirks['reaction_fp'] = smirks['reaction'].apply(
                lambda r: get_fp(r, fp=fp_type, concatenate=concatenate)
            )

            # plot similarity heatmap of all smirks
            smirks_similarities: pd.DataFrame = get_similarity_heatmap(smirks)
            smirks_similarities['fp_type'] = fp_type
            smirks_similarities['concatenate'] = concatenate
            smirks_sim_heatmaps.append(smirks_similarities)

            results_df, best_config = evaluate_w_hyperparameters(test_data, smirks, fp_type=fp_type,
                                                                 concatenate=concatenate)

            results_df['fp_type'] = fp_type
            results_df['concatenate'] = concatenate

            all_results.append(results_df)
    all_results_df = pd.concat(all_results)
    best_idx = all_results_df['accuracy'].idxmax()
    overall_best = all_results_df.iloc[best_idx]
    logger.info(f"Overall best configuration: {overall_best.to_dict()}")
    all_results_df.to_csv('../data/all_fingerprint_configurations.csv', index=False)
    path_to_fig = plot_overall_results(all_results_df)
    logger.info(f"Plot overall results: {path_to_fig}")
    paths = plot_heatmaps(smirks_sim_heatmaps)

    return overall_best


def main():
    test_data: pd.DataFrame = load_data(path='../data/2025-04-04_rxn_classifer_91_test_data_labeled.json')
    with open('../../syndirella/constants/RXN_SMIRKS_CONSTANTS.json', 'r') as f:
        data = json.load(f)
    smirks = pd.DataFrame(list(data.items()), columns=['name', 'reaction'])

    best_config = run_classifier(test_data=test_data, smirks=smirks)

    logger.info("\nClassifier Evaluation Results:")
    logger.info("-" * 50)
    logger.info(f"Best configuration:")
    logger.info(f"  Fingerprint type: {best_config['fp_type']}")
    logger.info(f"  Concatenate: {best_config['concatenate']}")
    logger.info(f"  Similarity metric: {best_config['similarity_metric']}")
    logger.info(f"  Top N: {best_config['top_n']}")
    logger.info(f"  Threshold: {best_config['threshold']}")
    logger.info(f"  Accuracy: {best_config['accuracy']:.4f}")


if __name__ == "__main__":
    main()
