# function splitSegment 
# splitSegment returns a list of segment boarder of given numerical attribute and classLable
# 
# Usage: 
#   splitSegment(attr,classLabel), returns segment boarder split points
# Arguments: 
#   attr: attrbute to be split
#   classLable: according class label of attribute, only for 2 class lable case
splitSegment <- function(attr,classLabel,minleaf = 1) {
    # sort the attribute list
    sortedAttr <- sort(unique(attr))
    # calculate all numeric split points
    lengthOfSortedAttr <- length(sortedAttr)
    if (lengthOfSortedAttr == 1) {
        return(NA)
    }
    splitPoints <- (sortedAttr[1 : lengthOfSortedAttr - 1] + sortedAttr[2 : lengthOfSortedAttr])/2
    if (minleaf != 1) {return(splitPoints)}
    # calculate the probability in class(1) for each item in the sortedAttr 
    prob <- numeric(lengthOfSortedAttr)
    for (i in 1:length(sortedAttr)) {
        probTemp <- classLabel[attr == sortedAttr[i]]
        prob[i] = sum(probTemp)/length(probTemp)
    }
    # the segment boarder split points are those points where the probablity changes
    segmentBoarder <- splitPoints[prob[1:length(prob) - 1] != prob[2:length(prob)]]
    return(segmentBoarder)
}

# function impurityGini
# impurity returns a impurity value of given node class label with Gini-index impurity function
impurityGini <- function(x) {
  length_x <- length(x)
  p1 <- sum(x)/length_x
  impurity <- p1 * (1 - p1)
  return(impurity)
}

# function resubstitutionErr
# return the resubstituotion error of a given node
# Arguments:
#   classLabel : the class label of samples in the node
resubstitutionError <- function(classLabel) {
    amountSamples <- length(classLabel)
    amountClass1 <- sum(classLabel)
    amountClass0 <- amountSamples - amountClass1
    resubstitutionError <- min(amountClass0,amountClass1)/amountSamples
    return(resubstitutionError)
} 

# function impurityReduction
# calculate the impurity reduction of a given split point on a given node
# arguments:
#   attr :      the attributes used to split the samples
#   splitPoint: the split point on the given attributes
#   classLabel: the class label of samples in the node
#   minleaf:    minmum amount of samples in a child node after splitting
impurityReduction <- function(attr,classLabel,splitPoint,minleaf = 1) {
    leftchildClassLabel <- classLabel[attr <= splitPoint]
    rightchildClassLabel <- classLabel[attr > splitPoint]
    # if the amount of samples in left or righ child node then set the impurity reduction to 0 which 
    # means it is a optional split point no matter what the real impurity reduction of the split is
    if ((length(leftchildClassLabel) < minleaf) || (length(rightchildClassLabel) < minleaf)) {
        impurityReduction <- 0
    }
    # calculate the impurity reduction with the formula delta(i) = i(t) - {p(l)i(l) + p(r)i(r)}
    else {
        impurityNode <- impurityGini(classLabel)
        impurityleftchildNode <- impurityGini(leftchildClassLabel)
        impurityrightchildNode <- impurityGini(rightchildClassLabel)
        impurityReduction <- impurityNode -
                             length(leftchildClassLabel)/length(classLabel) * impurityleftchildNode -
                             length(rightchildClassLabel)/length(classLabel) * impurityrightchildNode
    }
    return(impurityReduction)
}


# function bestSplit
# find out the best split points for a given node and all attributes and class label of the samples
# arguments:
#   data:       the full data set including all attributes and class label of all samples
#   sampleList: the list of the index of samples for current node
#   nmin:       the minimum amount of samples for a node to split
#   minleaf:    minimum amout of samples in a child node after splitting
bestSplit <- function(data,sampleList,nmin = 2,minleaf = 1) {
    # directly return a node with splitvar = -1(invalid) and splitval = -1.0(invalid)
    # if the number of samples in current node is less than nmin
    if (length(sampleList) < nmin) {
        node <- data.frame(parent = NA, lchild = NA, rchild = NA, splitvar = -1, splitval = -1.0, n = NA, imp = NA, gr = NA, pnode = NA, enode = NA, rnode = NA)
        return(node)
    }
    # get the class label for current node
    classLabel <- data[sampleList,ncol(data)]
    # create a vector of split point and according impurity reduction for each attribute
    splitVals <- numeric(ncol(data) - 1)
    deltaImpurityAllAttrs <- numeric(ncol(data) - 1)
    # loop for each attributes
    for (i in 1 : (ncol(data) - 1)) {
        attr <- data[sampleList,i]
        splitPoints <- splitSegment(attr,classLabel,minleaf)
        # if there is no split point in current attribute set impurity reduction to 0 and split point value to -1.0(invalid)
        if (is.na(splitPoints[1]) == TRUE) {
            deltaImpurityAllAttrs[i] = 0
            splitVals[i] = -1.0
        }
        else {
            deltaImpurity <- numeric(length(splitPoints))
            # loop for each split point in a same attribute
            for (j in 1: length(splitPoints)) {
                deltaImpurity[j] <- impurityReduction(attr,classLabel,splitPoints[j],minleaf)
            }
            # get the best split point in a single attribute
            maxDeltaImpurity <- max(deltaImpurity)
            deltaImpurityAllAttrs[i] <- maxDeltaImpurity
            splitVals[i] <- splitPoints[match(maxDeltaImpurity,deltaImpurity)]
        }
    }
    # get the best split point accross all attributes
    maxDeltaImpurityAllAttrs <- max(deltaImpurityAllAttrs)
    if (maxDeltaImpurityAllAttrs == 0) {
        splitVar <- -1
        splitVal <- -1.0
    }
    else {
        splitVal <- splitVals[match(maxDeltaImpurityAllAttrs,deltaImpurityAllAttrs)]
        splitVar <- match(splitVal,splitVals)
    }
    node <- data.frame(parent = NA, lchild = NA, rchild = NA, splitvar = splitVar, splitval = splitVal, n = NA, imp = NA, gr = NA, pnode = NA, enode = NA, rnode = NA, tieclass = NA)
    return(node)
}

# function tree.grow
# generate a classification tree with give samples with attributions and class labels
# arguments:
#   data:       the full data set including all attributes and class label of all samples
#   nmin:       the minimum amount of samples for a node to split
#   minleaf:    minimum amout of samples in a child node after splitting
tree.grow <- function(data,nmin = 2, minleaf = 1) {
    # create a empty tree
    tree <- data.frame()
    # create the sample list for root node
    sampleList <- 1:nrow(data)
    # pass the list of root node to node list
    # the format of a single load is a list of {samples and parent node No.}
    # So, the format of nodeList is a list of list
    nodeList <- list(list(samples = sampleList, parent = 0))
    # set the No. of root node to 1
    nodeNum <- 1
    # if the node list is not NULL, then grow the tree
    while (length(nodeList) != 0) {
        # get the current node from node list
        currNodeSampleList <- nodeList[[1]][['samples']]
        nodeParent <- nodeList[[1]][['parent']]
        # nodelist <- nodelist - current node
        nodeList[1] <- NULL
        # get the class label of all samples in current node
        currNodeSampleLabels <- data[currNodeSampleList,ncol(data)]
        # calculate the impurity for current node
        impurityCurrNode <- impurityGini(currNodeSampleLabels)
        # calculate the node class lable according marjority rule
        if (sum(currNodeSampleLabels) == length(currNodeSampleLabels)/2) {
            tieClass <- TRUE
            currNodeClass <- 0
        }
        else if (sum(currNodeSampleLabels) < length(currNodeSampleLabels)/2) {
            tieClass <- FALSE
            currNodeClass <- 0
            }
        else {
            tieClass <- FALSE
            currNodeClass <- 1
        }
        # get the length of current node
        nCurrNode <- length(currNodeSampleList)
        # calculate the best split of current node 
        node <- bestSplit(data,currNodeSampleList,nmin,minleaf)
        node['parent'] = nodeParent
        # if current node can be split, the apply the split
        if (impurityCurrNode > 0 && node[['splitvar']] != -1) {
            leftChildNodeSampleList <- currNodeSampleList[data[currNodeSampleList,node[['splitvar']]] <= node[['splitval']]]
            rightChildNodeSampleList <- currNodeSampleList[data[currNodeSampleList,node[['splitvar']]] > node[['splitval']]]
            nodeList <- append(nodeList,list(list(samples = leftChildNodeSampleList, parent = nodeNum)))
            node['lchild'] <- nodeNum + length(nodeList)
            nodeList <- append(nodeList,list(list(samples = rightChildNodeSampleList, parent = nodeNum)))
            node['rchild'] <- nodeNum + length(nodeList)
        }
        # else no clild node for current node
        else {
            node['lchild'] <- -1
            node['rchild'] <- -1
        }
        # calculate other parameters for current node
        node['n'] <- nCurrNode
        node['imp'] <- impurityCurrNode
        node['gr'] <- currNodeClass
        node['pnode'] <- length(currNodeSampleLabels) / nrow(data)
        node['enode'] <- resubstitutionError(currNodeSampleLabels)
        node['rnode'] <- node['pnode'] * node['enode']
        node['tieclass'] <- tieClass
        tree <- rbind(tree,node)
        nodeNum <- nodeNum + 1
    }
    return(tree)
}

# function tree.classify
# return the predicted class label of inputed samples with a given tree
# arguments:
#   attrs:  the attributes of inputed samples
#   tree:   inputed classification tree
tree.classify <- function(attrs,tree) {
    classLabels <- NULL
    # predict the class for each input sample
    for (i in 1:nrow(attrs)) {
        nodeNo = 1
        # search the tree according the input attribute value of the sample until leaf node reached
        while (nodeNo != -1) {
            leafNodeNo <- nodeNo
            splitVar = tree[[nodeNo,'splitvar']]
            splitVal = tree[[nodeNo,'splitval']]
            # set nodeNo to -1 to quit if the leaf node reached
            if (splitVar == -1) {
                nodeNo = -1
            }
            else {
                if (attrs[i,splitVar] <= splitVal) {
                    nodeNo <- tree[nodeNo,'lchild']
                }
                else {
                    nodeNo <- tree[nodeNo,'rchild']
                }
            }
        }
        # get the class lable for current sample and append it into the predicted lables
        classLabel <- tree[leafNodeNo,'gr']
        classLabels <- append(classLabels,classLabel)
    }
    return(classLabels)
}

# function tree.simplify
# simplify a classification tree with pruning the redudant leaf nodes untill there is no leaf node which has same
# majority class with sibling node or class tie in either or both siblings node
# arguments:
#   tree:   the tree object to be simplified
tree.simplify <- function(tree) {
    simplifiedTree <- tree
    canSimplify <- TRUE
    while (canSimplify) {
        canSimplify <- FALSE
        for (i in 1:nrow(tree)) {
            if (simplifiedTree[i,'parent'] != -1 && simplifiedTree[i,'splitvar'] == -1) {
                parentNode <- simplifiedTree[i,'parent']
                siblings <- c(simplifiedTree[parentNode,'lchild'],simplifiedTree[parentNode,'rchild'])
                sibling <- siblings[siblings != i]
                if (simplifiedTree[sibling,'parent'] != -1 && simplifiedTree[sibling,'splitvar'] == -1 && 
                    (simplifiedTree[i,'gr'] == simplifiedTree[sibling,'gr'] || 
                     simplifiedTree[i,'tieclass'] == TRUE || simplifiedTree[sibling,'tieclass'] == TRUE)) {
                    # simplify the tree
                    simplifiedTree[siblings,] = -1
                    simplifiedTree[parentNode,c('lchild','rchild','splitvar','splitval')] = -1
                    canSimplify <- TRUE
                }
            }
        }
    }
    return(simplifiedTree)
}


# function performanceTest
# return the confusion matrix of true class and predicted class
# arguments:
#   trueClass       : the true class label of samples
#   predictedClass  : the predicted class label of samples  
performanceTest <- function(trueClass,predictedClass) {
    if (length(trueClass) != length(predictedClass)) {
        return('ERROR! The length of inputed true class and predicted class are not equal!')
    }
    # t0p0 is the counter of the amount of samples with class label '0' and be predicted as '0', likewise for t0p1,....
    t0p0 <- 0
    t0p1 <- 0
    t1p0 <- 0
    t1p1 <- 0
    for (i in 1:length(trueClass)) {
        if (trueClass[i] == 0 && predictedClass[i] == 0) {
            t0p0 <- t0p0 + 1
        }
        else if (trueClass[i] == 0 && predictedClass[i] == 1) {
            t0p1 <- t0p1 + 1
        }
        else if (trueClass[i] == 1 && predictedClass[i] == 0) {
            t1p0 <- t1p0 + 1
        }
        else {
            t1p1 <- t1p1 + 1
        }
    }
    confusionMatrix <- matrix(c(t0p0,t0p1,t1p0,t1p1),2,2,byrow = TRUE,dimnames = list(c('True 0','True 1'),c('Predict 0','Predict 1')))
    return(confusionMatrix)
}

main <- function(dataFileName,header = TRUE, nmin = 2,minleaf = 1) {
    data <- readData(dataFileName,header)
    tree <- tree.grow(data,nmin, minleaf)
    simplifiedTree <- tree.simplify(tree)
    predictedClass <- tree.classify(data,simplifiedTree)
    confusionMatrix <- performanceTest(data[,ncol(data)],predictedClass)
    return(confusionMatrix)
}

train <- function(data) {
    heartDisease.dat <- data
    samples <- 1:nrow(heartDisease.dat)
    trainSamples <- sample(nrow(heartDisease.dat),200)
    testSamples <- samples[-trainSamples]
    
    trainSampleFolds <- list()
    for (i in 1:10) {
        fold <- sample(trainSamples,20)
        trainSampleFolds <- append(trainSampleFolds,list(fold))
        trainSamples <- trainSamples[-match(fold,trainSamples)]
    }
    
    nmin = 10
    minleaf = 3
    trainSamples <- samples[-testSamples]
    predictedClass <- NULL
    trueClass <- NULL
    for (i in 1:10) {
        predictSamples <- trainSampleFolds[[i]]
        growTreeSamples <- trainSamples[-match(predictSamples,trainSamples)]
        tree <- tree.grow(heartDisease.dat[growTreeSamples,],nmin,minleaf)
        simplifiedTree <- tree.simplify(tree)
        predictedClass <- append(predictedClass, tree.classify(heartDisease.dat[predictSamples,],simplifiedTree))
        trueClass <- append(trueClass, heartDisease.dat[predictSamples,ncol(heartDisease.dat)])
    }
    confusionMatrix <- performanceTest(trueClass,predictedClass)
    print(confusionMatrix)
}


