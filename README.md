# CliffordSim
 Efficient Clifford Simulator using the tableau method from Aaronson & Gottesman's
 [`Improved Simulation of Stabilizer Circuits`](https://arxiv.org/pdf/quant-ph/0406196) paper.

 Uses a base `Tableau` class to which gates are applied through methods such as `tab.CZ(c,t)`, `tab.H(i)`, `tab.S(i)`, etc.
