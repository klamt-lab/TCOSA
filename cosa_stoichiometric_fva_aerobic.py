from cosa_stoichiometric_fva import cosa_single_swap_test
import ray
ray.init(log_to_driver=False)

for base_nadx_scenario in ("WILDTYPE", "FLEXIBLE", "SINGLE_COFACTOR"):
    cosa_single_swap_test(
        anaerobic=False,
        reac_id="TEST_0_818",
        mu=0.818,
        base_nadx_scenario=base_nadx_scenario
    )
