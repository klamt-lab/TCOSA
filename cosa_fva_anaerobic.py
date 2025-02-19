"""Run of thermodynamics-using FVAs under anaerobic conditions and glucose."""

from cosa_fva import cosa_single_swap_test
import ray
ray.init(log_to_driver=False)


for base_nadx_scenario in ("WILDTYPE", "FLEXIBLE", "SINGLE_COFACTOR"):
    # cosa_single_swap_test(
    #     anaerobic=True,
    #     reac_id="TEST_0_371",
    #     mu=0.371,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    cosa_single_swap_test(
        anaerobic=True,
        reac_id="TEST_0_321",
        mu=0.321,
        base_nadx_scenario=base_nadx_scenario
    )
    # cosa_single_swap_test(
    #     anaerobic=True,
    #     reac_id="TEST_0_271",
    #     mu=0.271,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    # cosa_single_swap_test(
    #     anaerobic=True,
    #     reac_id="TEST_0_221",
    #     mu=0.221,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    # cosa_single_swap_test(
    #     anaerobic=True,
    #     reac_id="TEST_0_171",
    #     mu=0.171,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    # cosa_single_swap_test(
    #     anaerobic=True,
    #     reac_id="TEST_0_121",
    #     mu=0.121,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    # cosa_single_swap_test(
    #     anaerobic=True,
    #     reac_id="TEST_0_071",
    #     mu=0.071,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    # cosa_single_swap_test(
    #     anaerobic=True,
    #     reac_id="TEST_0_05",
    #     mu=0.05,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    # cosa_single_swap_test(
    #     anaerobic=True,
    #     reac_id="TEST_0_03",
    #     mu=0.05,
    #     base_nadx_scenario=base_nadx_scenario
    # )
