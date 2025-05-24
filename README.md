# Symbolic-Tensor-manipulation
This is a project that performs symbolically written operations between tensors. Some of these implemented operations are:  scalar contraction between indeces, standard matrix multiplication, or outer product between indeces.  This library was developed with relativity equation manipulation in mind.

# Syntax of symbolic operations

  1. Define tensor shape and initialize it with CreateTensor().
  2. initialize tensor components with constants with InitializeTensorConstant() or with function pointers with InitializeTensorFunction().
  3. Define symbolic expression of tensor equation. For example, if we want to represent the equation $A_{l} B_{jk} C_{kl} = D_{j}$, we do it writing the operation "l,jk,kl->j".
  4. Use the Einsum() function to perform the operation described symbolically with initialized tensors. 

# Usage

A few use cases are included in the main() function to prove library functionality:

  1. Parsing and identifying arbitrary symbolic index expressions with many inputs. Some operations under test: ijk,kl,li->j (operation: inner product between 3-tensor and matrix), ijk,kl->ijl (operation: inner product between 3-tensor and matrix), ij,k->ijk (operation: outer product between 2-tensor and vector), ij,ij-> (operation: scalar contraction between matrices).
  2. Numeric tensor operations (constant components). Different expressions are initialized numerically and the tensor result of the operation is printed in console using the Einsum() function.
  3. Numeric tensor initialization (variable components, using function pointers). Different tensor components are initialized with function pointers and are printed in console.
