# Fast Matrix Multiplication

This is a collection of snippets related to verifying small matrix multiplication formulas published in 2022, and provides conversions of them to the format used by THE FMM catalog at https://fmm.univ-lille.fr/

## DeepMind AlphaTensor
The directory [alphatensor](alphatensor) contains Maple scripts to verify the matrix
multiplication formulas discovered by the DeepMind team's AlphaTensor project. https://github.com/deepmind/alphatensor

## Kausers Moosbauer 2022

[KMtoFMM.mpl](KMtoFMM.mpl) and [KMtoFMM.ipynb](KMtoFFM.ipynb) are a script and a notebook that verify the
matrix multiplication formulas first reported by Kauers and Moosbauer in https://arxiv.org/abs/2210.04045

[KMtoFFM-AllFResults.mpl](KMtoFFM-AllFResults.mpl) is a script properly generalized to rectangular matrixes in arbitrary characteristic that verifies all the results reported in their Flip Graph paper https://arxiv.org/abs/2212.01175 and https://github.com/jakobmoosbauer/flips

