# python

import datetime
import json
import logging
import os

import glob2
import pandas as pd

# Import the class to test
from aizynth.AiZynthManager import AiZynthManager


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


def get_test_data(all_routes: dict, logger: logging.Logger) -> list[dict]:
    """Stores all test data from routes"""
    all_data = []
    logger.info(f'Writing reactions from {len(all_routes)} input molecules.')
    for smiles, manager in all_routes.items():
        if not hasattr(manager, 'routes'):
            all_reactions: list[dict] = [{'smiles': smiles}]  # no routes found
        else:
            all_reactions: list[dict] = manager.export_reactions_to_dict(manager.routes)
        all_data.extend(all_reactions)
    logger.info(f'Exported {len(all_data)} reactions.')
    return all_data


def export_test_data(all_test_data: list[dict], logger: logging.Logger, output_path_short: str = 'test_data') -> bool:
    """Exports test data to .json"""
    output_path: str = f'{output_path_short}_{len(all_test_data)}_test_data.json'
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(all_test_data, f, ensure_ascii=False, indent=4)
    logger.info(f'Wrote {len(all_test_data)} reactions to {output_path}')


def setup(logger: logging.Logger) -> AiZynthManager:
    logger.info(f"AIZYNTH config file set to {os.environ['AIZYNTH_CONFIG_FILE']}")
    manager = AiZynthManager()
    manager.logger = logger
    return manager


def sample_master(master: pd.DataFrame, n_smiles: int, logger: logging.Logger) -> pd.DataFrame:
    """Sample master with routes to test."""
    seeds_not_to_include = [0, 20]  # manual look over results
    for seed in range(0, 100):
        if seed not in seeds_not_to_include:
            sampled: pd.DataFrame = master.sample(n=n_smiles, random_state=seed)
            try:
                assert len(sampled[(sampled['1_reaction'].notna()) & sampled[
                    '2_reaction'].isna()]) >= 0.5 * n_smiles  # at least 50% are single step
                assert sampled['error_type'].str.contains('NoSynthesisRoute').any()
                assert len(sampled[(sampled['1_reaction'].notna()) & sampled['2_reaction'].notna()]) > 0
                logger.info(f'Sampled {len(sampled)} reactions with seed {seed}.')
                return sampled
            except AssertionError:
                logger.info(f'Sampling master df with seed {seed} did not generate required samples, retrying...')
                continue
    logger.error('Failed to sample master df with required conditions after 100 attempts.')
    return pd.DataFrame()  # Return an empty DataFrame if all attempts fail


def get_smiles(n_smiles: int, logger: logging.Logger, output_path_short: str) -> list[str]:
    """Gets a list of smiles to test from Syndirella outputs"""
    to_sample: list[str] = glob2.glob('../data/*.csv')
    logger.info(f'Found {len(to_sample)} output smiles to test.')
    all_dfs = []
    for path in to_sample:
        df = pd.read_csv(path)
        all_dfs.append(df)
    master: pd.DataFrame = pd.concat(all_dfs)
    master = master.dropna(subset=['smiles'])
    master = master[~(master['error_type'] == 'ScaffoldPlacementError')]  # remove all with placement error
    logger.info(f'Found {len(master)} potential routes to sample from.')
    columns_to_keep = [f'{num}_{suffix}' for num in range(1, 6) for suffix in ['r1_smiles', 'r2_smiles', 'reaction']]
    master = master[['smiles', 'inchi_key', 'error_type', 'error_message', 'route_uuid',
                     'hit1'] + columns_to_keep]
    sampled: pd.DataFrame = sample_master(master=master, n_smiles=n_smiles, logger=logger)
    output_path: str = f'{output_path_short}_{len(sampled)}_sampled.csv'
    sampled.to_csv(output_path, index=False)
    logger.info(f'Wrote {len(sampled)} potential routes to {output_path}')
    return list(sampled['smiles'])


def main():
    logger: logging.Logger = setup_logging()
    output_path_short: str = f'../data/{datetime.datetime.now().strftime("%Y-%m-%d")}_rxn_classifer'
    smiles: list[str] = get_smiles(n_smiles=5, logger=logger, output_path_short=output_path_short)
    logger.info(f'Running test data with {len(smiles)} smiles')
    all_routes: dict[str, AiZynthManager] = {}
    for smile in smiles:
        manager: AiZynthManager = setup(logger)
        manager.find_routes(target_smiles=smile)
        all_routes[smile]: AiZynthManager = manager
    all_test_data: dict = get_test_data(all_routes, logger)
    export_test_data(all_test_data=all_test_data, output_path_short=output_path_short, logger=logger)


if __name__ == "__main__":
    main()
