# DustNGwBBP

Study impact of dust NG on inferences on r

### How to run the pipeline

1. Update `config.yaml` with the parameters you want to use. Update the paths too.
2. Update `bbpower_input/write_sh_file.sh` with the same value for these parameters.
3. Make sure `utils/bbpw_script.py` also contains the setup you want to use for BBPower.
2. Run `test.py`.
3. Go to where BBPower is install in your system and run `write_sh_file.sh`.

This generates the final chains, which you can plot on your own. 

