# Each example must run and reproduce the tables and figures in examples/output
import os
import subprocess
import sys
import numpy as np
import pandas as pd
import pytest

root=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..")
expected=os.path.join(root,"examples","output")
# Example: decimals of its tables
examples={"Current":2,"Dragonera":2,"EstacaDeBares":2,"Reanalysis":3,
          "Synthetic":4,"WRF":2,"Ensembles":1}

@pytest.mark.parametrize("name",sorted(examples))
def test_example(name,tmp_path):
    subprocess.run([sys.executable,os.path.join(root,"examples",name+".py"),str(tmp_path)],
                   check=True,capture_output=True)
    produced=sorted(os.listdir(tmp_path))
    assert produced==sorted(f for f in os.listdir(expected) if f.startswith(name+"_"))
    for f in produced:
        if not f.endswith(".csv"):
            continue
        new=pd.read_csv(os.path.join(tmp_path,f))
        old=pd.read_csv(os.path.join(expected,f))
        assert list(new.columns)==list(old.columns)
        assert list(new["modelName"])==list(old["modelName"])
        # Allow one unit in the last decimal (rounding on other platforms)
        a=new.iloc[:,1:].to_numpy(dtype=float)
        b=old.iloc[:,1:].to_numpy(dtype=float)
        assert np.allclose(a,b,rtol=0.,atol=1.01*10.**-examples[name],equal_nan=True), f
