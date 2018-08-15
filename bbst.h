// C program to find number of elements
// greater than a given value in AVL
#include <stdio.h>
#include <stdlib.h>

struct TNode {
    float key;
    struct TNode* left, *right;
    int height;
    int desc;
    int count;
};
 
int height(struct TNode* N)
{
    if (N == NULL)
        return 0;
    return N->height;
}
 
// A utility function to get maximum
// of two integers
int max(int a, int b)
{
    return (a > b) ? a : b;
}
 
struct TNode* newNode(float key)
{
    struct TNode* node = (struct TNode*)
                    malloc(sizeof(struct TNode));
    node->key = key;
    node->left = NULL;
    node->right = NULL;
    node->height = 1; // initially added at leaf
    node->desc = 0;
    node->count = 1;
    return (node);
}
 
// A utility function to right rotate subtree
// rooted with y
struct TNode* rightRotate(struct TNode* y)
{
    struct TNode* x = y->left;
    struct TNode* T2 = x->right;
 
    // Perform rotation
    x->right = y;
    y->left = T2;
 
    // Update heights
    y->height = max(height(y->left), height(y->right)) + 1;
    x->height = max(height(x->left), height(x->right)) + 1;
 
    // calculate the number of children of x and y
    // which are changed due to rotation.
    int val = (T2 != NULL) ? T2->desc : -1;
    y->desc = y->desc - (x->desc + 1) + (val + 1);
    x->desc = x->desc - (val + 1) + (y->desc + 1);
 
    return x;
}
 
// A utility function to left rotate subtree rooted
// with x
struct TNode* leftRotate(struct TNode* x)
{
    struct TNode* y = x->right;
    struct TNode* T2 = y->left;
 
    // Perform rotation
    y->left = x;
    x->right = T2;
 
    // Update heights
    x->height = max(height(x->left), height(x->right)) + 1;
    y->height = max(height(y->left), height(y->right)) + 1;
 
    // calculate the number of children of x and y
    // which are changed due to rotation.
    int val = (T2 != NULL) ? T2->desc : -1;
    x->desc = x->desc - (y->desc + 1) + (val + 1);
    y->desc = y->desc - (val + 1) + (x->desc + 1);
 
    return y;
}
 
// Get Balance factor of node N
int getBalance(struct TNode* N)
{
    if (N == NULL)
        return 0;
    return height(N->left) - height(N->right);
}
 
struct TNode* insert(struct TNode* node, float key)
{
    /* 1. Perform the normal BST rotation */
    if (node == NULL)
        return (newNode(key));

    // If key already exists in BST, icnrement count and return
    if (key == node->key){
        node->count++;
        return node;
    }
 
    if (key < node->key) {
        node->left = insert(node->left, key);
        node->desc++;
    }
 
    else if (key > node->key) {
        node->right = insert(node->right, key);
        node->desc++;
    }
 
    /* 2. Update height of this ancestor node */
    node->height = 1 + max(height(node->left),
                           height(node->right));
 
    /* 3. Get the balance factor of this ancestor
        node to check whether this node became
        unbalanced */
    int balance = getBalance(node);
 
    // If node becomes unbalanced, 4 cases arise
 
    // Left Left Case
    if (balance > 1 && key < node->left->key)
        return rightRotate(node);
 
    // Right Right Case
    if (balance < -1 && key > node->right->key)
        return leftRotate(node);
 
    // Left Right Case
    if (balance > 1 && key > node->left->key) {
        node->left = leftRotate(node->left);
        return rightRotate(node);
    }
 
    // Right Left Case
    if (balance < -1 && key < node->right->key) {
        node->right = rightRotate(node->right);
        return leftRotate(node);
    }
 
    /* return the (unchanged) node pointer */
    return node;
}
 
/* Given a non-empty binary search tree, return the
   node with minimum key value found in that tree.
   Note that the entire tree does not need to be
   searched. */
struct TNode* minValueNode(struct TNode* node)
{
    struct TNode* current = node;
 
    /* loop down to find the leftmost leaf */
    while (current->left != NULL)
        current = current->left;
 
    return current;
}
 
// Recursive function to delete a node with given key
// from subtree with given root. It returns root of
// the modified subtree.
struct TNode* deleteNode(struct TNode* root, float key)
{
    // STEP 1: PERFORM STANDARD BST DELETE
 
    if (root == NULL)
        return root;
 
    // If the key to be deleted is smaller than the
    // root's key, then it lies in left subtree
    if (key < root->key) {
        root->left = deleteNode(root->left, key);
        root->desc = root->desc - 1;
    }
 
    // If the key to be deleted is greater than the
    // root's key, then it lies in right subtree
    else if (key > root->key) {
        root->right = deleteNode(root->right, key);
        root->desc = root->desc - 1;
    }
 
    // if key is same as root's key, then This is
    // the node to be deleted
    else {
        // If key is present more than once, simply decrement
        // count and return
        if (root->count > 1)
        {
            (root->count)--;
            return root;
        }
        // ElSE, delete the node
        // node with only one child or no child
        if ((root->left == NULL) || (root->right == NULL)) {
 
            struct TNode* temp = root->left ? 
                                root->left : root->right;
 
            // No child case
            if (temp == NULL) {
                temp = root;
                root = NULL;
                free(temp);
 
            } 
            else // One child case
            {
                *root = *temp; // Copy the contents of
                               // the non-empty child
                free(temp);
            }
        } else {
            // node with two children: Get the inorder
            // successor (smallest in the right subtree)
            struct TNode* temp = minValueNode(root->right);
 
            // Copy the inorder successor's data to this node
            root->key = temp->key;
 
            // Delete the inorder successor
            root->right = deleteNode(root->right, temp->key);
            root->desc = root->desc - 1;
        }
    }
 
    // If the tree had only one node then return
    if (root == NULL)
        return root;
 
    // STEP 2: UPDATE HEIGHT OF THE CURRENT NODE
    root->height = 1 + max(height(root->left),
                           height(root->right));
 
    // STEP 3: GET THE BALANCE FACTOR OF THIS NODE (to
    // check whether this node became unbalanced)
    int balance = getBalance(root);
 
    // If this node becomes unbalanced, 4 cases arise
 
    // Left Left Case
    if (balance > 1 && getBalance(root->left) >= 0)
        return rightRotate(root);
 
    // Left Right Case
    if (balance > 1 && getBalance(root->left) < 0) {
        root->left = leftRotate(root->left);
        return rightRotate(root);
    }
 
    // Right Right Case
    if (balance < -1 && getBalance(root->right) <= 0)
        return leftRotate(root);
 
    // Right Left Case
    if (balance < -1 && getBalance(root->right) > 0) {
        root->right = rightRotate(root->right);
        return leftRotate(root);
    }
 
    return root;
}

// A utility function to print preorder traversal of
// the tree.
void preOrder(struct TNode* root)
{
    if (root != NULL) {
        printf("%f (%d, d%d), ", root->key, root->count, root->desc);
        preOrder(root->left);
        preOrder(root->right);
    }
}
 
// // Returns count of
// int CountGreater(struct TNode* root, float x)
// {
//     int res = 0;
 
//     // Search for x. While searching, keep
//     // updating res if x is greater than
//     // current node.
//     while (root != NULL) {
 
//         int desc = (root->right != NULL) ? 
//                    root->right->desc : -1;
 
//         if (root->key > x) {
//             res = res + desc + 1 + 1;
//             root = root->left;
//         } else if (root->key < x)
//             root = root->right;
//         else {
//             res = res + desc + 1;
//             break;
//         }
//     }
//     return res;
// }

// Returns count of
int CountLesser(struct TNode* root, float x)
{
    int res = 0;

    if(root == NULL)
        return res;

    struct TNode* max_node=root; 
    while(max_node->right != NULL)
        max_node = max_node->right;

    if(x > max_node->key)
        return root->desc + root->count;
 
    // Search for x. While searching, keep
    // updating res if x is lesser than
    // current node.
    while (root != NULL) {
        int desc = (root->left != NULL) ? 
                   root->left->desc : -1;
        //printf("%d\n", root->count);
        if (root->key < x) {
            //printf("%d ", desc+2);
            res += desc + 1 + 1;
            root = root->right;
        } else if (root->key > x)
            root = root->left;
        else {
            res += desc + 1;
            break;
        }
    }
    return res;
} 

// /* Driver program to test above function*/
// int main()
// {
//     struct TNode* root = NULL;
//     root = insert(root, 0.0);
//     root = insert(root, 0.0);
//     root = insert(root, 2.0);
 
//     /* The constructed AVL Tree would be
//           9
//         /   \
//         1   10
//        / \     \
//        0  5     11
//       /  / \
//     -1  2   6    */
 
//     printf("Preorder traversal of the constructed AVL "
//            "tree is \n");
//     preOrder(root);
//     printf("\nNumber of elements lesser than 110.10 are %d",
//            CountLesser(root, 2.5));
 
//     root = deleteNode(root, 0.0);
 
//     /* The AVL Tree after deletion of 10
//          1
//         / \
//         0 9
//        /  / \
//      -1  5  11
//         / \
//         2  6 */
 
//     printf("\nPreorder traversal after deletion of 0.10 \n");
//     preOrder(root);
//     printf("\nNumber of elements lesser than 110.1000 are %d",
//            CountLesser(root, 5.5));
//     return 0;
// }
