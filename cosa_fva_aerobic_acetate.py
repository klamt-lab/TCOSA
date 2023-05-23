"""Run of thermodynamics-using FVAs under aerobic conditions and acetate."""

from cosa_fva import cosa_single_swap_test
import ray
ray.init(log_to_driver=False)

for base_nadx_scenario in ("WILDTYPE", "FLEXIBLE", "SINGLE_COFACTOR"):
    # cosa_single_swap_test(
    #     anaerobic=False,
    #     reac_id="TEST_0_206",
    #     mu=0.2059,
    #     base_nadx_scenario=base_nadx_scenario,
    #     c_source="acetate",
    # )
    cosa_single_swap_test(
        anaerobic=False,
        reac_id="TEST_0_181",
        mu=0.181,
        base_nadx_scenario=base_nadx_scenario,
        c_source="acetate",
    )
    # cosa_single_swap_test(
    #     anaerobic=False,
    #     reac_id="TEST_0_106",
    #     mu=0.106,
    #     base_nadx_scenario=base_nadx_scenario,
    #     c_source="acetate",
    # )
    # cosa_single_swap_test(
    #     anaerobic=False,
    #     reac_id="TEST_0_056",
    #     mu=0.056,
    #     base_nadx_scenario=base_nadx_scenario,
    #     c_source="acetate",
    # )
    # cosa_single_swap_test(
    #     anaerobic=False,
    #     reac_id="TEST_0_003",
    #     mu=0.03,
    #     base_nadx_scenario=base_nadx_scenario,
    #     c_source="acetate",
    # )
