#ifndef _EINSUM_H_
#define _EINSUM_H_

/*****************************************************************************************************************************
 *                                                          Includes
 *****************************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*****************************************************************************************************************************
 *                                                      Macro definitions
 *****************************************************************************************************************************/

#define SCALAR_CONTRACTION_1        0
#define SCALAR_CONTRACTION_2        1
#define OUTER_PRODUCT_1             2
#define MATRIX_MULTIPLICATION       3
#define VECTOR_CONTRACTION_1        4
#define OUTER_PRODUCT_2             5
#define OUTER_PRODUCT_3             6
#define OUTER_PRODUCT_4             7
#define INNER_PRODUCT_2             8
#define INNER_PRODUCT_3             9
#define INNER_PRODUCT_4             10
#define INNER_PRODUCT_5             11
#define INNER_PRODUCT_6             12


/*****************************************************************************************************************************
 *                                                   Expression declarations
 *****************************************************************************************************************************/

typedef struct {
    char *subscript;            // expression 
    size_t* input_shapes;       // object input shapes
    size_t num_inputs;          // number of inputs in expression
    size_t output_shape;        // shape of tensor object
    int operation_type;         // operation type to perform

} Expression;

/* ParseExpression
 * [in] : char*; expression to parse. Must be in format "(object 1 index),(object 2 indeces)->(output indeces)". Example: "ijk,jk->i"
 * [out] : Expression; expression object created after expression is parsed
 * 
 * Parses a string into an expression. Expression object holds number of input, operations to perform between tensors, among other things.
 */
Expression ParseExpression(char* notation);

/* FreeExpression
 * [in] : Expression*
 * [out] : void
 * 
 * Frees dynamic memory allocated for Expression* object.
 */
void FreeExpression(Expression* result);

/* PrintExpression
 * [in] : Expression*;
 * [out] : void
 * 
 * Prints Expression to console, reporting the value of all its atributes that were generated when operation was parsed.
 */
void PrintExpression(Expression* result);

/* GenerateOutputIndeces
 * [in] : char*; Expression to parse. Must be in format "(object 1 index),(object 2 indeces)->(output indeces)"
 * [out] : char*; Indeces of output object
 * 
 * Function parses inigial char* expression and automatically indicates what must be the output indeces of resulting object.
 */
char* GenerateOutputIndeces(char* expression);

/* Simplify
 * [in] : char*; Expression to parse.
 * [out] : Expression*
 * 
 * Parses an expression that can have many inputs (with one output) and returns an array of expression equations in the format 
 * "(object 1 index),(object 2 indeces)->(output indeces)". Example: For input "ij,jk,kl->il", the function would return the
 * following Expressions: "ij,jk->ik", "ik,kl->il".
 */
Expression* Simplify(char* expression);


/*****************************************************************************************************************************
 *                                                      Tensor declarations
 *****************************************************************************************************************************/

typedef struct {
    double* data;           // tensor data
    int dataSize;           // size (in components) of data array
    size_t* shape;          // shape of tensor
    size_t ndim;            // tensor dimension
    int isConstant;         // flag to indicate if a tensor is constant
} Tensor;

typedef double(*functionPointer)(double*);      // typedef definition for function pointers

/* CreateTensor
 * [in] : size_t, size_t; Tensor shape and tensor dimension
 * [out] : Tensor
 * 
 * Allocates dinamic memory necessary to hold all elements of a Tensor and returns said Tensor.
 */
Tensor CreateTensor(size_t* shape, size_t ndim);

/* InitializetensorConstant
 * [in] : Tensor*, double*; Tensor passed by reference, double array to initialize all its components with constant values
 * [out] : void
 * 
 * Initializes all (or some) of Tensor components with constant double values.
 */
void InitializeTensorConstant(Tensor *t, double* vals);

/* InitializeTensorFunction
 * [in] : Tensor*, functionPointer*, double*; Tensor passed by reference, function pointers to initialize each component, parameters to evaluate the function pointers
 * [out] : void
 * 
 * Initializes all (or some) of Tensor components with function values evaluated with certain parameters.
 */
void InitializeTensorFunction(Tensor *t, functionPointer* funcs, double* params);

/* PrintTensor
 * [in] : Tensor*; Tensor passed by reference
 * [out] : void
 * 
 * Prints information about the Tensor and prints all its components in console.
 */
void PrintTensor(const Tensor* t);

/* GetElement
 * [in] : Tensor*, size_t*; Tensor passed by reference, position of the desired value
 * [out] : double; Desired element
 * 
 * Gets the element located in the desired Tensor indeces.
 */
double GetElement(const Tensor* t, size_t* position);

/* SetElement
 * [in] : Tensor*, double, size_t*; Tensor passed by reference, value to set, position of value to set
 * [out] : void
 * 
 * Sets the element located in the desired Tensor indeces to a certain value.
 */
void SetElement(Tensor* t, double value, size_t* position);

/* Einsum
 * [in] : char*, Tensor*, size_t; operation to perform, Tensor array of inputs, number of operands in the Tensor array
 * [out] : Tensor
 * 
 * Uses the symbolic expression to calculate a result using the operands passed to the function and returns the Tensor result.
 */
Tensor Einsum(char* operation, Tensor* operands, size_t num_operands);

/* FreeTensor
 * [in] : Tensor*
 * [out] : void
 * 
 * Frees dynamic memory allocated for Tensor* object.
 */
void FreeTensor(Tensor* T);


/*****************************************************************************************************************************
 *                                                      Tensor Operations
 *****************************************************************************************************************************/

/* ScalarContraction1
 * [in] : Tensor*, Tensor*; Inputs passed by referece
 * [out] : Tensor
 * 
 * Calculates the following operations between inputs: "i,i->"
 */
Tensor ScalarContraction1(Tensor* A, Tensor* B);

/* ScalarContraction2
 * [in] : Tensor*, Tensor*; Inputs passed by referece
 * [out] : Tensor
 * 
 * Calculates the following operations between inputs: "ij,ij->"
 */
Tensor ScalarContraction2(Tensor* A, Tensor* B);

/* OuterProduct1
 * [in] : Tensor*, Tensor*; Inputs passed by referece
 * [out] : Tensor
 * 
 * Calculates the following operations between inputs: "i,j->ij"
 */
Tensor OuterProduct1(Tensor* A, Tensor* B);

/* MatrixMultiplications
 * [in] : Tensor*, Tensor*; Inputs passed by referece
 * [out] : Tensor
 * 
 * Calculates the following operations between inputs: "ij,jk->ik"
 */
Tensor MatrixMultiplication(Tensor* A, Tensor* B);

/* VectorContraction1
 * [in] : Tensor*, Tensor*; Inputs passed by referece
 * [out] : Tensor
 * 
 * Calculates the following operations between inputs: "ij,j->i"
 */
Tensor VectorContraction1(Tensor* A, Tensor* B); 

/* OuterProduct2
 * [in] : Tensor*, Tensor*; Inputs passed by referece
 * [out] : Tensor
 * 
 * Calculates the following operations between inputs: "ij,k->ijk"
 */
Tensor OuterProduct2(Tensor* A, Tensor* B);

/* OuterProduct3
 * [in] : Tensor*, Tensor*; Inputs passed by referece
 * [out] : Tensor
 * 
 * Calculates the following operations between inputs: "ij,kl->ijkl"
 */
Tensor OuterProduct3(Tensor* A, Tensor* B);

/* OuterProduct4
 * [in] : Tensor*, Tensor*; Inputs passed by referece
 * [out] : Tensor
 * 
 * Calculates the following operations between inputs: "ijk,l->ijkl"
 */
Tensor OuterProduct4(Tensor* A, Tensor* B);

/* InnerProduct2
 * [in] : Tensor*, Tensor*; Inputs passed by referece
 * [out] : Tensor
 * 
 * Calculates the following operations between inputs: "ijk,k->ij"
 */
Tensor InnerProduct2(Tensor* A, Tensor* B);

/* InnerProduct3
 * [in] : Tensor*, Tensor*; Inputs passed by referece
 * [out] : Tensor
 * 
 * Calculates the following operations between inputs: "ijk,kl->ijl"
 */
Tensor InnerProduct3(Tensor* A, Tensor* B);

/* InnerProduct4
 * [in] : Tensor*, Tensor*; Inputs passed by referece
 * [out] : Tensor
 * 
 * Calculates the following operations between inputs: "ijkl,kl->ij"
 */
Tensor InnerProduct4(Tensor* A, Tensor* B);

/* InnerProduct5
 * [in] : Tensor*, Tensor*; Inputs passed by referece
 * [out] : Tensor
 * 
 * Calculates the following operations between inputs: "ijk,jk->i"
 */
Tensor InnerProduct5(Tensor* A, Tensor* B);

/* InnerProduct5
 * [in] : Tensor*, Tensor*; Inputs passed by referece
 * [out] : Tensor
 * 
 * Calculates the following operations between inputs: "ijkl,l->ijk"
 */
Tensor InnerProduct6(Tensor* A, Tensor* B);


/*****************************************************************************************************************************
 *                                                        Filesaving
 *****************************************************************************************************************************/

/* SaveTensorEquation
 * [in] : char*, Tensor*, functionPointer*, int, size_t*, size_t, double*, int, double, int, char*; Operation to perform,
 * Tensor operands, array of functionPointers to initialize all operands, number of operands, output shape, output dimensions,
 * parameters, position of sweeped parameter, parameter final value, number of sweep points, savefile path
 * [out] : void
 * 
 * Calculates operations on tensor operands initialized with function pointers. A set parameter is sweeped and the result
 * of the tensor operation is saved in a .txt file. Format for the text file is: Each row contains all resulting Tensor elements
 * for a given sweeped parameter value and each column is just a given Tensor element for different parameter sweep values
 */
void SaveTensorEquation(char* expression, Tensor* operands, functionPointer *funcs, int numOperands, size_t *output_shape, 
                        size_t output_ndim, double* params, int paramPosition, double paramEnd, int N, char* savePath);


#endif
