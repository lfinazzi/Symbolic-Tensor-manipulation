#include "Einsum.h"

Expression ParseExpression(char* notation) {

    Expression result;
    result.subscript = notation;

    // Find the arrow indicating output shape
    const char* arrow_pos = strstr(notation, "->");

    if (arrow_pos == NULL || arrow_pos == notation) {
        fprintf(stderr, "Error: Invalid einsum notation. Missing or misplaced '->'.\n");
        exit(EXIT_FAILURE);
    }

    // Extract input and output parts of the notation
    size_t input_part_length = arrow_pos - notation;
    size_t output_part_length = strlen(notation) - (arrow_pos - notation) - 2;

    char* input_part = (char*)malloc((input_part_length + 1) * sizeof(char));
    char* output_part = (char*)malloc((output_part_length + 1) * sizeof(char));

    strncpy(input_part, notation, input_part_length);
    input_part[input_part_length] = '\0';       // terminates the array with a string terminator

    strncpy(output_part, arrow_pos + 2, output_part_length);
    output_part[output_part_length] = '\0';       // terminates the array with a string terminator

    // Count the number of inputs
    size_t num_inputs = 0;
    char *input_part_copy = input_part;
    int i;
    for (i = 0; input_part_copy[i]; input_part_copy[i]==',' ? i++ : *input_part_copy++);

    if(i == 0){
        fprintf(stderr, "Error: Invalid einsum notation. Missing or misplaced ','.\n");
        exit(EXIT_FAILURE);
    }
    num_inputs = i + 1;

    // Extract input shapes
    size_t* input_shapes = (size_t*)malloc(num_inputs * sizeof(size_t));
    for(int i = 0; i < num_inputs; i++)
        input_shapes[i] = 0;

    int l = 0;
    for (const char* c = input_part; *c != '\0'; ++c) {
        if (*c >= 'a' && *c <= 'z') {
            input_shapes[l] += 1;
        }
        else if (*c == ',')
            l++;
    }

    // Extract output shape
    size_t num_output_indices = 0;
    for (const char* c = output_part; *c != '\0'; ++c) {
        if (*c >= 'a' && *c <= 'z') {
            num_output_indices++;
        }
    }

    // save in Expression object
    result.output_shape = num_output_indices;
    result.input_shapes = input_shapes;
    result.num_inputs = num_inputs;

    // Determine which operation should the Expression hold
    if(num_inputs == 2){

        if(num_output_indices == 0){
            if(result.input_shapes[0] == 1){
                result.operation_type = SCALAR_CONTRACTION_1;
            }
            else if(result.input_shapes[0] == 2){
                result.operation_type = SCALAR_CONTRACTION_2;
            }
        }

        else if(num_output_indices == 1){
            
            if((result.input_shapes[0] == 2 && result.input_shapes[1] == 1) ||
                        (result.input_shapes[0] == 1 && result.input_shapes[1] == 2)){
                result.operation_type = VECTOR_CONTRACTION_1;
            }
            else if((result.input_shapes[0] == 3 && result.input_shapes[1] == 2) ||
                        (result.input_shapes[0] == 2 && result.input_shapes[1] == 3)){
                result.operation_type = INNER_PRODUCT_5;
            }
        }

        else if(num_output_indices == 2){
            
            if(result.input_shapes[0] == 1 && result.input_shapes[1] == 1){
                result.operation_type = OUTER_PRODUCT_1;
            }
            else if(result.input_shapes[0] == 2){
                result.operation_type = MATRIX_MULTIPLICATION;
            }
            else if((result.input_shapes[0] == 3 && result.input_shapes[1] == 1) ||
                        (result.input_shapes[0] == 1 && result.input_shapes[1] == 3)){
                result.operation_type = INNER_PRODUCT_2;
            }
            else if((result.input_shapes[0] == 4 && result.input_shapes[1] == 2) ||
                        (result.input_shapes[0] == 2 && result.input_shapes[1] == 4)){
                result.operation_type = INNER_PRODUCT_4;
            }

        }
        
        else if(num_output_indices == 3){
            
            if((result.input_shapes[0] == 2 && result.input_shapes[1] == 1) ||
                        (result.input_shapes[0] == 1 && result.input_shapes[1] == 2)){
                result.operation_type = OUTER_PRODUCT_2;
            }           
            else if((result.input_shapes[0] == 3 && result.input_shapes[1] == 2) ||
                        (result.input_shapes[0] == 3 && result.input_shapes[1] == 2)){
                result.operation_type = INNER_PRODUCT_3;
            }
            else if((result.input_shapes[0] == 4 && result.input_shapes[1] == 1) ||
                        (result.input_shapes[0] == 1 && result.input_shapes[1] == 4)){
                result.operation_type = INNER_PRODUCT_6;
            }

        }

        else if(num_output_indices == 4){
            
            if((result.input_shapes[0] == 2 && result.input_shapes[1] == 2)){
                result.operation_type = OUTER_PRODUCT_3;
            }
            else if((result.input_shapes[0] == 3 && result.input_shapes[1] == 1) ||
                        (result.input_shapes[0] == 1 && result.input_shapes[1] == 3)){
                result.operation_type = OUTER_PRODUCT_4;
            }

        }

        else{       // operation assign error
            fprintf(stderr, "Couldn't assign an operation type.\n");
            exit(EXIT_FAILURE);
        }

    }

    return result;
}

void FreeExpression(Expression* result) {
    free(result->subscript);
    free(result->input_shapes);
}

void PrintExpression(Expression* result)
{
    // Information parsed in expression
    printf("\nExpression: %s\n", result->subscript);
    printf("Num inputs: %ld\n", result->num_inputs);
    printf("Input indeces: (");
    for(int i = 0; i < result->num_inputs - 1; i++)
        printf("%ld,", result->input_shapes[i]);
    printf("%ld", result->input_shapes[result->num_inputs - 1]);
    printf(")\n");
    printf("Output indeces: %ld\n", result->output_shape);

    if(result->num_inputs > 2){     // can't assign name to composite expression
        printf("Expression type: Composite expression with %d inputs\n\n", result->num_inputs);
        return;
    }

    // Expression type
    printf("Expression type: ");
    if(result->operation_type == SCALAR_CONTRACTION_1)
        printf("Scalar contraction between vectors.\n");
    else if(result->operation_type == SCALAR_CONTRACTION_2)
        printf("Scalar contraction between matrices.\n");
    else if(result->operation_type == OUTER_PRODUCT_1)
        printf("Outer product between vectors.\n");
    else if(result->operation_type == MATRIX_MULTIPLICATION)
        printf("Matrix multiplication.\n");
    else if(result->operation_type == VECTOR_CONTRACTION_1)
        printf("Vector contraction.\n");
    else if(result->operation_type == OUTER_PRODUCT_2)
        printf("Outer product between matrix and vector.\n");
    else if(result->operation_type == OUTER_PRODUCT_3)
        printf("Outer product between matrix and matrix.\n");
    else if(result->operation_type == OUTER_PRODUCT_4)
        printf("Outer product between 3-tensor and vector.\n");
    else if(result->operation_type == INNER_PRODUCT_2)
        printf("Outer product between 3-tensor and vector.\n");
    else if(result->operation_type == INNER_PRODUCT_3)
        printf("Inner product between 3-tensor and matrix.\n");
    else if(result->operation_type == INNER_PRODUCT_4)
        printf("Inner product between 4-tensor and matrix.\n");
    else if(result->operation_type == INNER_PRODUCT_5)
        printf("Inner product between 3-tensor and matrix.\n");
    else if(result->operation_type == INNER_PRODUCT_6)
        printf("Inner product between 4-tensor and vector.\n");
   
    printf("\n");
    return;
}

char* GenerateOutputIndeces(char* expression)
{
    char* out = (char*)malloc(5*sizeof(char));              // 5 is enough for any implemented output

    int i = 0;
    int k = 0;
    for (const char* c = expression; *c != '\0'; ++c) {     // iterate over all chars in expression
        char currChar = *c;
        int repeated = 0;

        for (const char* d = expression; *d != '\0'; ++d) { // see if character repeats twice (because it counts char with itself)
            if(*d == currChar)
                repeated += 1;
        }

        if(repeated < 2 && currChar != ',' && currChar != '-' && currChar != '>'){  // identify repeated char and generate output indeces
            out[k] = currChar;
            k++;
        }

        i++;
    }

    return out;
}

Expression* Simplify(char* expression)
{
    Expression original = ParseExpression(expression);

    // allocates all subscripts for all expressions
    char** inputs = (char**)malloc((original.num_inputs)*sizeof(char*));
    for(int i = 0; i < original.num_inputs; i++){
        inputs[i] = (char*)malloc(20*sizeof(char));      // 20 chars will always be enough for each operation
    }

    const char* arrow_pos = strstr(expression, "->");

    if (arrow_pos == NULL || arrow_pos == expression) {
        fprintf(stderr, "Error: Invalid einsum notation. Missing or misplaced '->'.\n");
        exit(EXIT_FAILURE);
    }

    // Extract input and output parts of the notation
    size_t input_part_length = arrow_pos - expression;
    char* input_part = (char*)malloc((input_part_length + 1) * sizeof(char));
    strncpy(input_part, expression, input_part_length);
    input_part[input_part_length] = '\0';

    int pos = 0;
    for (int i = 0; i < original.num_inputs; i++){
        strncpy(inputs[i], input_part + pos, original.input_shapes[i]);
        inputs[i][20] = '\0';
        pos += original.input_shapes[i] + 1;
    }

    // allocate memory for subscripts for each Expression
    char **simplifiedExpressions = (char**)malloc((original.num_inputs - 1)*sizeof(char*));
    for(int i = 0; i < original.num_inputs - 1; i++){
        simplifiedExpressions[i] = (char*)malloc(20*sizeof(char));      // 20 chars will always be enough for each operation
    }

    // Generate first operation subscript
    char *outputIndeces = (char*)malloc(20 * sizeof(char));

    strcat(simplifiedExpressions[0], inputs[0]);
    strcat(simplifiedExpressions[0], ",");
    strcat(simplifiedExpressions[0], inputs[1]);
    strcat(simplifiedExpressions[0], "->");

    outputIndeces = GenerateOutputIndeces(simplifiedExpressions[0]);
    strcat(simplifiedExpressions[0], outputIndeces);

    // Generate all the following subscript expressions using the previous subscript output
    // indeces as the first object
    for(int i = 1; i < original.num_inputs - 1; i++)
    {
        strcat(simplifiedExpressions[i], outputIndeces);
        strcat(simplifiedExpressions[i], ",");
        strcat(simplifiedExpressions[i], inputs[i + 1]);
        strcat(simplifiedExpressions[i], "->");

        outputIndeces = GenerateOutputIndeces(simplifiedExpressions[i]);
        strcat(simplifiedExpressions[i], outputIndeces);
    }
    
    Expression *exp = (Expression*)malloc((original.num_inputs - 1)*sizeof(Expression));

    // Create all output expressions
    for(int i = 0; i < original.num_inputs - 1; i++){
        exp[i] = ParseExpression(simplifiedExpressions[i]);
    }

    return exp;
}


/*****************************************************************************************************************************
 *                                                      Tensor definitions
 *****************************************************************************************************************************/

// Function to create a tensor with given shape
Tensor CreateTensor(size_t* shape, size_t ndim) {
    Tensor t;
    t.ndim = ndim;
    t.shape = (size_t*)malloc(ndim * sizeof(size_t));
    
    size_t size = 1;
    for (size_t i = 0; i < ndim; ++i) {
        t.shape[i] = shape[i];
        size *= shape[i];
    }

    t.data = (double*)malloc(size * sizeof(double));
    t.dataSize = size;
    t.isConstant = 0;

    return t;
}

void InitializeTensorConstant(Tensor *t, double* vals)
{
    t->data = vals; 
    return;
}

void InitializeTensorFunction(Tensor *t, functionPointer* funcs, double* params)
{
    int N = t->dataSize;
    for(int i = 0; i < N; i++)
        t->data[i] = funcs[i](params);
    t->isConstant = 0;      // false
    return;
}

// Function to print the tensor in a matrix-like format
void PrintTensor(const Tensor *t) {
    
    if(t->ndim == 0){
        printf("Matrix-like print of tensor with shape (%zu):\n", t->shape[0]);
        printf("%f ", t->data[0]);
        printf("\n");
    }    
    else if(t->ndim == 1){
        printf("Matrix-like print of tensor with shape (%zu):\n", t->shape[0]);

        for (size_t i = 0; i < t->shape[0]; ++i) {
                printf("%f ", t->data[i]);
            }
            printf("\n");
    }
    else if(t->ndim == 2){
        printf("Matrix-like print of tensor with shape (%zu, %zu):\n", t->shape[0], t->shape[1]);

        for (size_t i = 0; i < t->shape[0]; ++i) {
            for (size_t j = 0; j < t->shape[1]; ++j) {
                printf("%f ", GetElement(t, (size_t[]){i, j}));
            }
            printf("\n");
        }
        printf("\n");
    }
    else if(t->ndim == 3){
        printf("Matrix-like print of tensor with shape (%zu, %zu, %zu):\n", t->shape[0], t->shape[1], t->shape[2]);

        for (size_t i = 0; i < t->shape[0]; ++i) {
            for (size_t j = 0; j < t->shape[1]; ++j) {
                for (size_t k = 0; k < t->shape[2]; ++k) {
                    printf("%f ", GetElement(t, (size_t[]){i, j, k}));
                }
                printf("\n");
            }
            printf("\n");
        }
    }
    else if(t->ndim == 4){
        printf("Matrix-like print of tensor with shape (%zu, %zu, %zu, %zu):\n", t->shape[0], t->shape[1], t->shape[2], t->shape[3]);

        for (size_t i = 0; i < t->shape[0]; ++i) {
            for (size_t j = 0; j < t->shape[1]; ++j) {
                for (size_t k = 0; k < t->shape[2]; ++k) {
                    for (size_t l = 0; l < t->shape[3]; ++l) {
                        printf("%f ", GetElement(t, (size_t[]){i, j, k, l}));
                    }
                    printf("\n");
                }
                printf("\n");
            }
            printf("\n");
        }
    }
    else{
        fprintf(stderr, "Invalid Tensor dimension!\n");
        exit(EXIT_FAILURE);
    }
}

// Function to access elements like a matrix
double GetElement(const Tensor* t, size_t* pos) {
    if(t->ndim == 1)
        return t->data[pos[0]];
    else if(t->ndim == 2)
        return t->data[pos[0] * t->shape[1] + pos[1]];
    else if(t->ndim == 3)
        return t->data[pos[0] * (t->shape[1] * t->shape[2]) + pos[1] * t->shape[2] + pos[2]];
    else if(t->ndim == 4)
        return t->data[pos[0] * (t->shape[1] * t->shape[2] * t->shape[3]) + pos[1] * (t->shape[2] * t->shape[3]) + pos[2] * (t->shape[3]) + pos[3]];
    else{
        fprintf(stderr, "Can't fetch Tensor element\n");
        exit(EXIT_FAILURE);
    }
}

// Function to access elements like a matrix
void SetElement(Tensor* t, double value, size_t* pos) {
    if(t->ndim == 1){
        t->data[pos[0]] = value;
        return;
    }
    else if(t->ndim == 2){
        t->data[pos[0] * t->shape[1] + pos[1]] = value;
        return;
    }
    else if(t->ndim == 3){
        t->data[pos[0] * (t->shape[1] * t->shape[2]) + pos[1] * t->shape[2] + pos[2]] = value;
        return;
    }    
    else if(t->ndim == 4){
        t->data[pos[0] * (t->shape[1] * t->shape[2] * t->shape[3]) + pos[1] * (t->shape[2] * t->shape[3]) + pos[2] * (t->shape[3]) + pos[3]] = value;
        return;
    }
    fprintf(stderr, "Can't set Tensor element\n");
    exit(EXIT_FAILURE);
}

Tensor Einsum(char* operation, Tensor* operands, size_t num_operands) {

    Expression *exps = Simplify(operation);     // Generate all Expressions to evaluate sequentially

    Expression expression;

    Tensor result;
    Tensor resultCalc;
    for(int i = 0; i < num_operands - 1; i++)
    {
        expression = exps[i];
        if(i == 0){
            resultCalc = operands[i];
            result = resultCalc;
        }
        resultCalc = result;
        
        if(expression.operation_type == SCALAR_CONTRACTION_1){                                // vector dot product
            result = ScalarContraction1(&resultCalc, &operands[i + 1]);
        }
        else if(expression.operation_type == SCALAR_CONTRACTION_2){                           // matrix dot product
            result = ScalarContraction2(&resultCalc, &operands[i + 1]);
        }        
        else if(expression.operation_type == OUTER_PRODUCT_1){                                // vector outer product
            result = OuterProduct1(&resultCalc, &operands[i + 1]);
        }
        else if(expression.operation_type == MATRIX_MULTIPLICATION){                          // matrix multiplication
            result = MatrixMultiplication(&resultCalc, &operands[i + 1]);
        }
        else if(expression.operation_type == VECTOR_CONTRACTION_1){                           // vector contraction
            result = VectorContraction1(&resultCalc, &operands[i + 1]);
        }
        else if(expression.operation_type == OUTER_PRODUCT_2){                                // tensor outer product 1
            result = OuterProduct2(&resultCalc, &operands[i + 1]);
        }
        else if(expression.operation_type == OUTER_PRODUCT_3){                                // tensor outer product 2
            result = OuterProduct3(&resultCalc, &operands[i + 1]);
        }
        else if(expression.operation_type == OUTER_PRODUCT_4){                                // tensor outer product 3
            result = OuterProduct4(&resultCalc, &operands[i + 1]);
        }
        else if(expression.operation_type == INNER_PRODUCT_2){                                // tensor inner product 2
            result = InnerProduct2(&resultCalc, &operands[i + 1]);
        }
        else if(expression.operation_type == INNER_PRODUCT_3){                                // tensor inner product 3
            result = InnerProduct3(&resultCalc, &operands[i + 1]);
        }
        else if(expression.operation_type == INNER_PRODUCT_4){                               // tensor inner product 4
            result = InnerProduct4(&resultCalc, &operands[i + 1]);
        }
        else if(expression.operation_type == INNER_PRODUCT_5){                               // tensor inner product 5
            result = InnerProduct5(&resultCalc, &operands[i + 1]);
        }
        else if(expression.operation_type == INNER_PRODUCT_6){                               // tensor inner product 5
            result = InnerProduct6(&resultCalc, &operands[i + 1]);
        }
        else{
            fprintf(stderr, "Can't assign operation!\n");
            exit(EXIT_FAILURE);
        }
    }

    return result;
}

void FreeTensor(Tensor* T)
{
    free(T->shape);
    free(T->data);
    return;
}

/*****************************************************************************************************************************
 *                                                 Tensor Operation definitions
 *****************************************************************************************************************************/

// i,i->
Tensor ScalarContraction1(Tensor* A, Tensor* B)
{
    Tensor result = CreateTensor((size_t[]){1}, 0);           // scalar

    result.data[0] = 0;
    for(size_t i = 0; i < A->shape[0]; i++){
        result.data[0] += A->data[i] * B->data[i];
    }
    return result;
}

// ij,ij->
Tensor ScalarContraction2(Tensor* A, Tensor* B)
{
    Tensor result = CreateTensor((size_t[]){1}, 0);           // scalar
    Tensor mult = MatrixMultiplication(A, B);
    
    result.data[0] = 0;
    for(size_t i = 0; i < mult.shape[0]; i++){
        result.data[0] += GetElement(&mult, (size_t[]){i, i});
    }
    return result;
}

// i,j->ij
Tensor OuterProduct1(Tensor* A, Tensor* B)
{
    Tensor result = CreateTensor((size_t[]){A->shape[0], B->shape[0]}, 2);      // 2-tensor

    for (size_t i = 0; i < result.shape[0]; ++i) {
        for (size_t j = 0; j < result.shape[1]; ++j) {
            SetElement(&result, A->data[i] * B->data[j], (size_t[]){i, j}); 
        }
    }
    return result;
}

// ij,jk->ik
Tensor MatrixMultiplication(Tensor* A, Tensor* B)
{
    Tensor result = CreateTensor((size_t[]){A->shape[0], B->shape[1]}, 2);      // 2-tensor
    double value;

    // Perform the contraction (sum over common index)
    for (size_t i = 0; i < result.shape[0]; ++i) {
        for (size_t k = 0; k < result.shape[1]; ++k) {
            value = 0;
            for (size_t j = 0; j < A->shape[1]; ++j) {
                value += GetElement(A, (size_t[]){i, j}) * GetElement(B, (size_t[]){j, k});
            }
            SetElement(&result, value, (size_t[]){i, k}); 
        }
    }
    return result;
    
}

// ij,j->i
Tensor VectorContraction1(Tensor* A, Tensor* B)
{
    Tensor result = CreateTensor((size_t[]){A->shape[0]}, 1);                   // vector

    // Perform the contraction (sum over common index)
    for (size_t i = 0; i < result.shape[0]; ++i) {
            result.data[i] = 0;
            for (size_t k = 0; k < B->shape[0]; ++k) {
                result.data[i] += GetElement(A, (size_t[]){i, k}) * B->data[k];
            }
    }
    return result;    

}

// ij,k->ijk
Tensor OuterProduct2(Tensor* A, Tensor* B)
{
    Tensor result = CreateTensor((size_t[]){A->shape[0], A->shape[1], B->shape[0]}, 3);          // 3-tensor
    double value = 0;

    for (size_t i = 0; i < result.shape[0]; ++i) {
        for (size_t j = 0; j < result.shape[1]; ++j) {
            for (size_t k = 0; k < result.shape[2]; ++k) {
                value = GetElement(A, (size_t[]){i, j}) * B->data[k];
                SetElement(&result, value, (size_t[]){i, j, k});
            }
        }
    }
    return result;    

}

// ij,kl->ijkl
Tensor OuterProduct3(Tensor* A, Tensor* B)
{
    Tensor result = CreateTensor((size_t[]){A->shape[0], A->shape[1], B->shape[0], B->shape[1]}, 4);          // 4-tensor
    double value = 0;

    for (size_t i = 0; i < result.shape[0]; ++i) {
        for (size_t j = 0; j < result.shape[1]; ++j) {
            for (size_t k = 0; k < result.shape[2]; ++k) {
                for (size_t l = 0; l < result.shape[3]; ++l) {
                    value = GetElement(A, (size_t[]){i, j}) * GetElement(B, (size_t[]){k, l});
                    SetElement(&result, value, (size_t[]){i, j, k, l});
                }
            }
        }
    }
    return result; 
}

// ijk,l->ijkl
Tensor OuterProduct4(Tensor* A, Tensor* B)
{
    Tensor result = CreateTensor((size_t[]){A->shape[0], A->shape[1], A->shape[2], B->shape[0]}, 4);          // 4-tensor
    double value = 0;

    for (size_t i = 0; i < result.shape[0]; ++i) {
        for (size_t j = 0; j < result.shape[1]; ++j) {
            for (size_t k = 0; k < result.shape[2]; ++k) {
                for (size_t l = 0; l < result.shape[3]; ++l) {
                    value = GetElement(A, (size_t[]){i, j, k}) * B->data[l];
                    SetElement(&result, value, (size_t[]){i, j, k, l});
                }
            }
        }
    }
    return result;    

}

// ijk,k->ij
Tensor InnerProduct2(Tensor* A, Tensor* B)
{
    Tensor result = CreateTensor((size_t[]){A->shape[0], A->shape[1]}, 2);          // matrix

    for (size_t i = 0; i < A->shape[0]; ++i) {
        for (size_t j = 0; j < A->shape[1]; ++j) {
            double value = 0;
            for (size_t k = 0; k < A->shape[2]; ++k) {
                value += GetElement(A, (size_t[]){i, j, k}) * B->data[k];
            }
            SetElement(&result, value, (size_t[]){i, j});
        }
    }
    return result; 
}

// ijk,kl->ijl
Tensor InnerProduct3(Tensor* A, Tensor* B)
{
    Tensor result = CreateTensor((size_t[]){A->shape[0], A->shape[1], B->shape[1]}, 3);          // 3-tensor
    double value = 0;

    for (size_t i = 0; i < A->shape[0]; ++i) {
        for (size_t j = 0; j < A->shape[1]; ++j) {
            for (size_t l = 0; l < B->shape[1]; ++l) {
                value = 0;
                for (size_t k = 0; k < A->shape[2]; ++k) {
                    value += GetElement(A, (size_t[]){i, j, k}) * GetElement(B, (size_t[]){k, l});
                }
                SetElement(&result, value, (size_t[]){i, j, l});
            }
        }
    }
    return result; 
}

// ijkl,kl->ij
Tensor InnerProduct4(Tensor* A, Tensor* B)
{
    Tensor result = CreateTensor((size_t[]){A->shape[0], A->shape[1]}, 2);          // matrix
    double value = 0;

    for (size_t i = 0; i < A->shape[0]; ++i) {
        for (size_t j = 0; j < A->shape[1]; ++j) {
            value = 0;
            for (size_t k = 0; k < A->shape[2]; ++k) {
                for (size_t l = 0; l < A->shape[3]; ++l) {
                    value += GetElement(A, (size_t[]){i, j, k, l}) * GetElement(B, (size_t[]){k, l});
                }  
            }
            SetElement(&result, value, (size_t[]){i, j});
        }
    }
    return result; 
}

// ijk,jk->i
Tensor InnerProduct5(Tensor* A, Tensor* B)
{
    Tensor result = CreateTensor((size_t[]){A->shape[0]}, 1);          // vector
    double value = 0;

    for (size_t i = 0; i < A->shape[0]; ++i) {
        value = 0;
        for (size_t j = 0; j < A->shape[1]; ++j) {
            for (size_t k = 0; k < A->shape[2]; ++k) {
                    value += GetElement(A, (size_t[]){i, j, k}) * GetElement(B, (size_t[]){j, k});
                }  
            }
            SetElement(&result, value, (size_t[]){i});
        }
    
    return result; 
}

// ijkl,l->ijk
Tensor InnerProduct6(Tensor* A, Tensor* B)
{
    Tensor result = CreateTensor((size_t[]){A->shape[0], A->shape[1], A->shape[2]}, 3);          // 3-tensor
    double value = 0;

    for (size_t i = 0; i < A->shape[0]; ++i) {
        for (size_t j = 0; j < A->shape[1]; ++j) {
            for (size_t k = 0; k < A->shape[2]; ++k) {
                value = 0;
                for (size_t l = 0; l < A->shape[3]; ++l) {
                    value += GetElement(A, (size_t[]){i, j, k, l}) * GetElement(B, (size_t[]){l});
                }
            SetElement(&result, value, (size_t[]){i, j, k});
            }
        }
    }
    return result; 
}
