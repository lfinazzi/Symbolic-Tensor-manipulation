#include "Einsum.h"
#include "functions.h"

int main(void)
{

    printf("*****************************\n");
    printf("Test program - Einsum library\n");
    printf("*****************************\n\n");

    // 1 ----------------------------------------------------------------
    // ------------------------------------------------------------------
    printf("1. parsing arbitrary expressions with many inputs:\n");
    char *operation = "ijk,kl,li->j";
    int numInputs = 3;

    printf("1.1:\n");
    printf("Original expression: %s\n", operation);
    
    Expression *ts = Simplify(operation);
    for(int i = 0; i < numInputs - 1; i++)
        PrintExpression(&ts[i]);
    
    operation = "ijkl,l,k->ij";

    printf("1.2:\n");
    printf("Original expression: %s\n", operation);
    
    ts = Simplify(operation);
    for(int i = 0; i < numInputs - 1; i++)
        PrintExpression(&ts[i]);

    operation = "ij,k,l,kl,ij->";
    numInputs = 5;
    printf("1.3:\n");
    printf("Original expression: %s\n", operation);
    
    ts = Simplify(operation);
    for(int i = 0; i < numInputs - 1; i++)
        PrintExpression(&ts[i]);
    
    // 2 ----------------------------------------------------------------
    // ------------------------------------------------------------------
    printf("2. Tensor operations, initialization with constants:\n\n");

    // Tensor definitions

    // square matrices for ease of evaluation
    size_t shape_a[] = {2, 2, 2, 2};
    size_t shape_b[] = {2, 2, 2};
    size_t shape_c[] = {2, 2};
    size_t shape_d[] = {2};
    // Allocates all tensor memory
    Tensor A = CreateTensor(shape_a, 4);  
    Tensor B = CreateTensor(shape_b, 3);  
    Tensor C = CreateTensor(shape_c, 2);  
    Tensor D = CreateTensor(shape_d, 1);  

    // initialize all Tensors with constant values
    InitializeTensorConstant(&A, (double[]){2, 3, 
                                            0, 1, 
                                            
                                            0, 2, 
                                            3, 1, 

                                            2, 1, 
                                            0, 0, 

                                            5, 2, 
                                            1, 1});

    InitializeTensorConstant(&B, (double[]){2, 3, 
                                            9, 1,

                                            0, 2,
                                            3, 1});

    InitializeTensorConstant(&C, (double[]){0, 1, 
                                            2, 3});

    InitializeTensorConstant(&D, (double[]){2, 3});

    operation = "ij,kl->ijkl";
    Expression expression = ParseExpression(operation);
    printf("2.1:\n");
    PrintExpression(&expression);
    Tensor E = Einsum(operation, (Tensor[]){C, C}, 2);
    PrintTensor(&E);

    operation = "ijkl,l->ijk";
    expression = ParseExpression(operation);
    printf("2.2:\n");
    PrintExpression(&expression);
    E = Einsum(operation, (Tensor[]){A, D}, 2);
    PrintTensor(&E);

    operation = "ij,jk->ik";
    expression = ParseExpression(operation);
    printf("2.3:\n");
    PrintExpression(&expression);
    E = Einsum(operation, (Tensor[]){C, C}, 2);
    PrintTensor(&E);

    operation = "ij,k,l->ijkl";
    expression = ParseExpression(operation);
    printf("2.4:\n");
    PrintExpression(&expression);
    E = Einsum(operation, (Tensor[]){C, D, D}, 3);
    PrintTensor(&E);

    operation = "ijk,k->ij";
    expression = ParseExpression(operation);
    printf("2.5:\n");
    PrintExpression(&expression);
    E = Einsum(operation, (Tensor[]){B, C}, 2);
    PrintTensor(&E);    

    operation = "ij,i,j->";
    expression = ParseExpression(operation);
    printf("2.6:\n");
    PrintExpression(&expression);
    E = Einsum(operation, (Tensor[]){C, D, D}, 3);
    PrintTensor(&E);   

    printf("\n");

    // 3 ----------------------------------------------------------------
    // ------------------------------------------------------------------
    printf("3. Tensor initialization with function pointers:\n\n");

    // Tensor definitions

    // square matrices for ease of evaluation
    size_t shape_f[] = {3, 3};
    // Allocates all tensor memory
    Tensor F = CreateTensor(shape_f, 2);

    // defines all functions to initialize allocated tensors
    functionPointer funcs[] = {fd, fnd1, fnd2,
                               fnd2, fd, fnd1,
                               fnd1, fnd2, fd}; 

    double params[] = {4, 3};

    printf("3.1:\n");
    InitializeTensorFunction(&F, funcs, params);
    PrintTensor(&F); 

    printf("3.2:\n");
    params[0] = 7;
    params[1] = 2;
    InitializeTensorFunction(&F, funcs, params);
    PrintTensor(&F); 

    return 0;
}