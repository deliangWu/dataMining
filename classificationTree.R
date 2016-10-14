# function tree.grow
# generate a classification tree (data.frame format) with give samples with attributions and class labels
# arguments:
#   data:       data.frame format, the full data set including all attributes and class label of all samples
#   nmin:       the minimum amount of samples for a node to split
#   minleaf:    minimum amout of samples in a child node after splitting
tree.grow <- function(x,y,nmin = 2, minleaf = 1) {
    # create a empty tree
    tree <- data.frame()
    # create the sample list for root node
    data <- cbind(x,y)
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
            } else if (sum(currNodeSampleLabels) < length(currNodeSampleLabels)/2) {
            tieClass <- FALSE
            currNodeClass <- 0
            } else {
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
        }else { # else no clild node for current node
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

# function tree.simplify
# simplify a classification tree with pruning the redundant leaf nodes until there is no leaf node which has same
# majority class with sibling node or class tie in either or both siblings node
# arguments:
#   tree:   data.format format, the tree object to be simplified
tree.simplify  <- function(tree){
    for (i in nrow(tree):1){   #from bottom
        if (tree[i,"lchild"]==-1) next   #skip leaf nodes
        if (tree[tree[i,"lchild"],"lchild"]==-1 && tree[tree[i,"rchild"],"lchild"]==-1 && tree[tree[i,"lchild"],"gr"]==tree[tree[i,"rchild"],"gr"]){
            tree[c(tree[i,"lchild"],tree[i,"rchild"]),]=-1
            tree[i,c('lchild','rchild','splitvar','splitval')]=-1
        }
    }
    tree<-removeUnusedNodes(tree)
    return(tree)
}

# function tree.classify
# return the predicted class label of inputed samples with a given tree
# arguments:
#   attrs:  data.frame format, the attributes of inputed samples
#   tree:   data.frame format, inputed classification tree
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
            }else {
                if (attrs[i,splitVar] <= splitVal) {
                    nodeNo <- tree[nodeNo,'lchild']
                }else {
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

# function tree.prune
# prune classifiction tree with Cost Complexity Pruning Algorithm, ruturn the pruned tree(data.frame format)
# arguments:
#   tree    : data.frame format, input classification tree before pruning
#   alpha   : penalty coefficient
tree.prune  <- function(tree,alpha=0){
    sampleList<-1:nrow(tree)
    for (i in nrow(tree):1){   #from bottom
        if (tree[i,"lchild"]==-1) next   #skip leaf nodes
        leafList<-as.vector(searchLeaf(i,tree,rep()))
        #calculate difference of complexity
        R <-( sum(tree[leafList,"enode"]*tree[leafList,"n"]) - tree[i,"enode"]*tree[i,"n"] )/ tree[1,"n"] + alpha*(length(leafList)-1)
        if (R >=0) {
            nodeList<-as.vector(searchNode(i,tree,rep()))
            nodeList<-nodeList[-1]
            tree[nodeList,]=-1
            tree[i,c('lchild','rchild','splitvar','splitval')]=-1
        }
    }
    tree<-removeUnusedNodes(tree)
    return(tree)
}


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

impurityEntropy <- function(x) {
  p1 <- sum(x)/length(x)
  impurity <- -p1 * log(p1) - (1 - p1) * log(1 - p1)
  return(impurity)
}

# function impurityGini
# impurity returns a impurity value of given node class label with Gini-index impurity function
impurityGini <- function(x) {
    #  length_x <- length(x)
    #  p1 <- sum(x)/length_x
    p1 <- sum(x)/length(x)
    impurity <- p1 * (1 - p1)
    #impurity <- -p1 * log(p1) - (1 - p1) * log(1 - p1)
    return(impurity)
}

# function resubstitutionErr
# return the resubstituotion error of a given node
# Arguments:
#   classLabel : the class label of samples in the node
resubstitutionError <- function(classLabel) {
    resubstitutionError <- min((length(classLabel) - sum(classLabel)),sum(classLabel))/length(classLabel)
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
    }else {# calculate the impurity reduction with the formula delta(i) = i(t) - {p(l)i(l) + p(r)i(r)}
        impurityNode <- impurityGini(classLabel)
        impurityleftchildNode <- impurityGini(leftchildClassLabel)
        impurityrightchildNode <- impurityGini(rightchildClassLabel)
        impurityReduction <- impurityNode -
                             length(leftchildClassLabel)/length(classLabel) * impurityleftchildNode -
                             length(rightchildClassLabel)/length(classLabel) * impurityrightchildNode
    }
    return(impurityReduction)
}


# function findDS
# used to replace loop structure in bestSplit function
# get the best split points in attributes (deltaImpurityAllAttr and splitVal)
# input: attribute(vector), ClassLabel(vector), minleaf(integer)
# output: a vector of deltaImpurityAllAttr and splitVal(vector)
findDS<-function(attr,classLabel,minleaf){
    #attr <- data[sampleList,i]
    splitPoints <- splitSegment(attr,classLabel,minleaf)
    # if there is no split point in current attribute set impurity reduction to 0 and split point value to -1.0(invalid)
    if (is.na(splitPoints[1]) == TRUE) {
        deltaImpurityAllAttr = 0
        splitVal = -1.0
    }else {
        # for each split point in a same attribute
        deltaImpurity<-apply(as.matrix(splitPoints),1,impurityReduction,"attr"=attr,"classLabel"=classLabel,"minleaf"=minleaf)
        
        # get the best split point in a single attribute
        maxDeltaImpurity <- max(deltaImpurity)
        deltaImpurityAllAttr <- maxDeltaImpurity
        splitVal<- splitPoints[match(maxDeltaImpurity,deltaImpurity)]
    }
    return(rbind(deltaImpurityAllAttr,splitVal))
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
    attr<-data[sampleList,1 : (ncol(data) - 1)]
    # loop for each attribute(use apply function instead)
    # to get the best split points in attributes
    DS<- apply(attr,2,findDS,"classLabel"=classLabel,"minleaf"=minleaf)
    deltaImpurityAllAttrs<-as.vector(DS[1,])
    splitVals<-as.vector(DS[2,])
    
    # get the best split point accross all attributes
    maxDeltaImpurityAllAttrs <- max(deltaImpurityAllAttrs)
    if (maxDeltaImpurityAllAttrs == 0) {
        splitVar <- -1
        splitVal <- -1.0
    }else {
        splitVal <- splitVals[match(maxDeltaImpurityAllAttrs,deltaImpurityAllAttrs)]
        splitVar <- match(splitVal,splitVals)
    }
    node <- data.frame(parent = NA, lchild = NA, rchild = NA, splitvar = splitVar, splitval = splitVal, n = NA, imp = NA, gr = NA, pnode = NA, enode = NA, rnode = NA, tieclass = NA)
    return(node)
}


searchLeaf<-function(nodeNr,tree,leafList){
    if (tree[nodeNr,"lchild"]==-1) {
        leafList<-cbind(leafList,nodeNr)
        return(leafList)
    }
    leafList<-searchLeaf(tree[nodeNr,"lchild"],tree,leafList)
    leafList<-searchLeaf(tree[nodeNr,"rchild"],tree,leafList)
    return(leafList)
}


searchNode<-function(nodeNr,tree,nodeList){
    nodeList<-cbind(nodeList,nodeNr)
    if (tree[nodeNr,"lchild"]==-1) {
        return(nodeList)
    }
    nodeList<-searchNode(tree[nodeNr,"lchild"],tree,nodeList)
    nodeList<-searchNode(tree[nodeNr,"rchild"],tree,nodeList)
    return(nodeList)
}


map <- function(x,sampleList) {
    if (!(x%in%sampleList)) return(x)
    which(sampleList==x)
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
#     # t0p0 is the counter of the amount of samples with class label '0' and be predicted as '0', likewise for t0p1,....
#     t0p0<-sum((1-trueClass)*(1-predictedClass))
#     t0p1<-sum((1-trueClass)*predictedClass)
#     t1p0<-sum(trueClass*(1-predictedClass))
#     t1p1<-sum(trueClass*predictedClass)
#     
#     confusionMatrix <- matrix(c(t0p0,t0p1,t1p0,t1p1),2,2,byrow = TRUE,dimnames = list(c('True 0','True 1'),c('Predict 0','Predict 1')))
    
    confusionMatrix<-tapply(rep(1,length(predictedClass)),list(trueClass,predictedClass),sum)
    confusionMatrix[is.na(confusionMatrix)]<-0
    
    return(confusionMatrix)
}

# no pruning，10次求平均
train <- function() {
    hddata <- read.csv("heartDisease.txt")
    samples <- 1:nrow(hddata)
    
    trainSamples <- sample(nrow(hddata),200)
    trainSamples0 <- trainSamples
    testSamples <- samples[-trainSamples]
    accu_array<-array(NA,cbind(80,40,1))
    
    for (k in 1:1){ #nr of trials
        trainSampleFolds <- list()
        trainSamples <- trainSamples0
        for (j in 1:10) {
            fold <- sample(trainSamples,20)
            trainSampleFolds <- append(trainSampleFolds,list(fold))
            trainSamples <- trainSamples[-match(fold,trainSamples)]
        }
        for (nmin in 2:80){
            for (minleaf in 1:min((nmin-1),40)){
                trainSamples <- trainSamples0
                predictedClass <- NULL
                trueClass <- NULL
                for (i in 1:10) {
                    predictSamples <- trainSampleFolds[[i]]
                    growTreeSamples <- trainSamples[-match(predictSamples,trainSamples)]
                    tree <- tree.grow(hddata[growTreeSamples,],nmin,minleaf)
                    simplifiedTree <- tree.simplify(tree)
                    predictedClass <- append(predictedClass, tree.classify(hddata[predictSamples,],simplifiedTree))
                    trueClass <- append(trueClass, hddata[predictSamples,ncol(hddata)])
                }
                confusionMatrix <- performanceTest(trueClass,predictedClass)
                accu_array[nmin,minleaf,k]<-(confusionMatrix[1,1]+confusionMatrix[2,2])/sum(confusionMatrix)
            }
        }
    }
    accu_array_kmean_new80<-accu_array
    accu<-apply(accu_array,c(1,2),mean)
    accu[is.na(accu)]<-0
    
    bestParameters<-which(accu==max(accu),arr.ind=T)
    bestNmin=as.numeric(bestParameters[1,1])
    bestMinleaf=as.numeric(bestParameters[1,2])
    
    bestTree <- tree.grow(hddata[trainSamples0,],bestNmin, bestMinleaf)
    bestSimplifiedTree <- tree.simplify(bestTree)
    bestPredictedClass <- tree.classify(hddata[testSamples,],bestSimplifiedTree)
    bestConfusionMatrix <- performanceTest(hddata[testSamples,ncol(hddata)],bestPredictedClass)
    bestAccu<-(bestConfusionMatrix[1,1]+bestConfusionMatrix[2,2])/sum(bestConfusionMatrix)
    bestAccu
    
    ######################################################
    
    
    accu<-1-apply(accu_array_kmean_new2,c(1,2),mean)
    
    # accu<-1-accu_array_kmean_new2[,,4]
    rts<-min(accu,na.rm=TRUE)
    se<-sqrt(rts*(1-rts)/200)
    rts+se
    which(accu<(rts+se),arr.ind = TRUE)
    which(accu==min(accu,na.rm=TRUE),arr.ind=TRUE)
    
    ##################################
    
    bestNmin=83
    bestMinleaf=28
    
    bestTree <- tree.grow(hddata[trainSamples0,],bestNmin, bestMinleaf)
    bestSimplifiedTree <- tree.simplify(bestTree)
    nrow(bestSimplifiedTree)
    bestPredictedClass <- tree.classify(hddata[testSamples,],bestSimplifiedTree)
    bestConfusionMatrix <- performanceTest(hddata[testSamples,ncol(hddata)],bestPredictedClass)
    bestAccu<-(bestConfusionMatrix[1,1]+bestConfusionMatrix[2,2])/sum(bestConfusionMatrix)
    1-bestAccu
    
    plotTree(bestSimplifiedTree,hddata)
}

train2 <- function(alpha=0) {
    
    accu_array<-array(0,cbind(2,1,101))
    for (j in 1:10) {
        fold <- sample(trainSamples,20)
        trainSampleFolds <- append(trainSampleFolds,list(fold))
        trainSamples <- trainSamples[-match(fold,trainSamples)]
    }
    for (k in 1:101){
        alpha= (k-1)*0.001
        for (nmin in 2:2){
            for (minleaf in 1:1){
                trainSamples <- trainSamples0
                predictedClass <- NULL
                trueClass <- NULL
                for (i in 1:10) {
                    predictSamples <- trainSampleFolds[[i]]
                    growTreeSamples <- trainSamples[-match(predictSamples,trainSamples)]
                    tree <- tree.grow(hddata[growTreeSamples,],nmin,minleaf)
                    simplifiedTree <- tree.simplify(tree)
                    prunedTree <- tree.prune (simplifiedTree, alpha)
                    predictedClass <- append(predictedClass, tree.classify(hddata[predictSamples,],prunedTree))
                    trueClass <- append(trueClass, hddata[predictSamples,ncol(hddata)])
                }
                confusionMatrix <- performanceTest(trueClass,predictedClass)
                accu_array[nmin,minleaf,k]<-(confusionMatrix[1,1]+confusionMatrix[2,2])/sum(confusionMatrix)
            }
        }
    }
    
    accu_alpha<-as.vector(accu_array[2,1,])
    
    bestParameters<-which(accu_array==max(accu_array),arr.ind=T)
    
    bestNmin=as.numeric(bestParameters[1,1])
    bestMinleaf=as.numeric(bestParameters[1,2])
    bestAlpha=as.numeric(bestParameters[1,3]-1)*0.001
    
    bestTree <- tree.grow(hddata[trainSamples0,],bestNmin, bestMinleaf)
    bestSimplifiedTree <- tree.simplify(bestTree)
    bestPrunedTree <- tree.prune (bestSimplifiedTree, bestAlpha)
    bestPredictedClass <- tree.classify(hddata[testSamples,],bestPrunedTree)
    bestConfusionMatrix <- performanceTest(hddata[testSamples,ncol(hddata)],bestPredictedClass)
    bestAccu<-(bestConfusionMatrix[1,1]+bestConfusionMatrix[2,2])/sum(bestConfusionMatrix)
    bestAccu
    
    #####################
    
    x<-1:101
    x<-0.001*(x-1)
    plot(x,accu_alpha,xlab="alpha",ylab="error rate",main="",pch=19)
    abline(h=rts,lty=3)
    abline(h=rts+se,lty=3)
    ############################
    accu_alpha<-(1-accu_alpha_new2)
    
    rts<-min(1-accu_alpha_new2)
    se<-sqrt(rts*(1-rts)/200)
    rts+se
    (which(accu_alpha<(rts+se))-1)*0.001
    (which(accu_alpha==min(accu_alpha,na.rm=TRUE))-1)*0.001
    ################################
    bestNmin=2
    bestMinleaf=1
    bestAlpha=0.011
    
    
    bestTree <- tree.grow(hddata[trainSamples0,],bestNmin, bestMinleaf)
    bestSimplifiedTree <- tree.simplify(bestTree)
    bestPrunedTree <- tree.prune (bestSimplifiedTree, bestAlpha)
    bestPredictedClass <- tree.classify(hddata[testSamples,],bestPrunedTree)
    bestConfusionMatrix <- performanceTest(hddata[testSamples,ncol(hddata)],bestPredictedClass)
    bestAccu<-(bestConfusionMatrix[1,1]+bestConfusionMatrix[2,2])/sum(bestConfusionMatrix)
    bestAccu
    
    plotTree(bestPrunedTree,hddata)
}

#bagging
train3 <- function(alpha=0) {
    
    bestAccu<-rep(0,10)
    for (i in 1:10){
        
        treeList <- bagging(hddata[trainSamples0,],1000,40)
        bestPredictedClass <- bagging.classify(hddata[testSamples,],treeList)
        bestConfusionMatrix <- performanceTest(hddata[testSamples,ncol
                                                      (hddata)],bestPredictedClass)
        bestAccu[i]<-(bestConfusionMatrix[1,1]+bestConfusionMatrix[2,2])/sum(bestConfusionMatrix)
    }
    
    bestAccu
    mean(1-bestAccu)
    sd(1-bestAccu)
}

#加入pruning
train1 <- function(alpha=0) {
    hddata <- read.csv("heartDisease.txt")
    samples <- 1:nrow(hddata)
    trainSamples <- sample(nrow(hddata),200)
    trainSamples0 <- trainSamples
    testSamples <- samples[-trainSamples]
    
    trainSampleFolds <- list()
    
    accu_array<-array(0,cbind(20,20,10))
    for (j in 1:10) {
        fold <- sample(trainSamples,20)
        trainSampleFolds <- append(trainSampleFolds,list(fold))
        trainSamples <- trainSamples[-match(fold,trainSamples)]
    }
    for (k in 1:10){
        alpha= (k-1)*0.005
        for (nmin in 1:20){
            for (minleaf in 1:20){
                trainSamples <- trainSamples0
                predictedClass <- NULL
                trueClass <- NULL
                for (i in 1:10) {
                    predictSamples <- trainSampleFolds[[i]]
                    growTreeSamples <- trainSamples[-match(predictSamples,trainSamples)]
                    tree <- tree.grow(hddata[growTreeSamples,],nmin,minleaf)
                    simplifiedTree <- tree.simplify(tree)
                    prunedTree <- tree.prune (simplifiedTree, alpha)
                    predictedClass <- append(predictedClass, tree.classify(hddata
                                                                           [predictSamples,],prunedTree))
                    trueClass <- append(trueClass, hddata[predictSamples,ncol(hddata)])
                }
                confusionMatrix <- performanceTest(trueClass,predictedClass)
                accu_array[nmin,minleaf,k]<-(confusionMatrix[1,1]+confusionMatrix[2,2])/sum
                (confusionMatrix)
            }
        }
    }
    
    bestParameters<-which(accu_array==max(accu_array),arr.ind=T)
    bestNmin=as.numeric(bestParameters[1,1])
    bestMinleaf=as.numeric(bestParameters[1,2])
    bestAlpha=as.numeric(bestParameters[1,3]-1)*0.005
    
    bestTree <- tree.grow(hddata[trainSamples0,],bestNmin, bestMinleaf)
    bestSimplifiedTree <- tree.simplify(bestTree)
    bestPrunedTree <- tree.prune (bestSimplifiedTree, bestAlpha)
    bestPredictedClass <- tree.classify(hddata[testSamples,],bestPrunedTree)
    bestConfusionMatrix <- performanceTest(hddata[testSamples,ncol
                                                  (hddata)],bestPredictedClass)
    bestAccu<-(bestConfusionMatrix[1,1]+bestConfusionMatrix[2,2])/sum(bestConfusionMatrix)
    bestAccu
}

#function removeUnusedNodes
#input: tree(dataframe)
#output: tree(dataframe) 
#remove the unused nodes (simplified or pruned, "parent"==-1), and reassign node nr. according to new row nr.
removeUnusedNodes<-function(tree){
    
    sampleList<-1:nrow(tree)
    sampleList<-sampleList[tree[,"parent"]!=-1]
    
    tree<-tree[!(tree[,"parent"]==-1),] #remove nodes
    
    #reassign node nr. in "parent","lchild","rchild", according to new row nr.
    for (i in cbind("parent","lchild","rchild")){
        treeNode<-as.matrix(tree[,i])
        treeNode<-apply(treeNode,1,map,"sampleList"=sampleList)
        tree[,i]<-treeNode
    }
    
    rownames(tree)<-1:nrow(tree)#reorder rownames
    return(tree)
}

#function bagging
#
#input:  data (data frame) : training samples
#        m (integer)       : how many trees would be constructed
#        n (integer)       : how many samples would be used to construct each tree
#output: treeList (list)   : a list of m trees
bagging<-function(data,m=100,n=25){
    data0<-data
    treeList<-list()
    sampleMatrix<-apply(as.matrix(1:m),1,sample,"x"=nrow(data0),"size"=n) #construct a sample matrix of n*m size
    for (i in 1:m){
        sample<-as.vector(sampleMatrix[,i])
        data<-data0[sample,]
        tree <- tree.grow(data,2, 1)
        simplifiedTree <- tree.simplify(tree)
        treeList[[i]]<-simplifiedTree
    }
    return(treeList)
}

#function bagging.classify
#
#input:  attrs (data frame)      : test samples
#        treeList (list)         : a list of m trees
#output: predictedClass (vector) : the predicted classes of test samples

bagging.classify<-function(attrs,treeList){
    predictedClassMatrix<-matrix(NA,length(treeList),nrow(attrs))
    
    for (i in 1:length(treeList)){
        predictedClassMatrix[i,]<-tree.classify(attrs,treeList[[i]])
    }
    predictedClass<-round(apply(predictedClassMatrix,2,mean))
    
    confusionMatrix <- performanceTest(attrs[,ncol(attrs)],predictedClass)
    accu<-(confusionMatrix[1,1]+confusionMatrix[2,2])/sum(confusionMatrix)
    accu
    
    return(predictedClass)

}
# function plotTree
# plot the input tree
# arguments:
#   tree            : inputed tree
#   data            : data used to generate the tree
plotTree <- function(tree,data,dispNodeInf = TRUE) {
    # load igraph library, need to install the igraph library
    library(igraph) 
    # get the name of attributes of the database
    attrNames <- names(data)
    # get all nodes except the root node
    node <- 2:nrow(tree)
    # get node's parent node
    parent <- tree$parent[-1]
    # generate edge lable for each edge
    parent.lchild <- tree[parent,'lchild']
    parent.splitvar <- tree[parent,'splitvar']
    parent.splitvar.name <- attrNames[parent.splitvar]
    parent.splitval <- tree[parent,'splitval']
    mux <- function(x) {if (x) {y <- ' < '} else {y <- ' >= '}}
    text <- apply(as.array(node == parent.lchild),1,mux)
    text <- paste0(parent.splitvar.name,text,parent.splitval)
    
    # generate tree graph
    dat <- data.frame(parent=parent,  node=node,  text=text)
    g <- graph.data.frame(dat)
    print(V(g))
    
    # get the node information including 'n', 'gr', 'enode' which will be used to display
    n <- tree[as.numeric(V(g)$name),'n']
    classLabel <- tree[as.numeric(V(g)$name),'gr']
    err <- round(tree[as.numeric(V(g)$name),'enode'],2)
    if (dispNodeInf)
        V(g)$name <- paste0("n=",n,' , err%=',err)
    # set the format of vertex
    V(g)$color <- ifelse(classLabel, 'pink','green')
    V(g)$label.cex <- 0.8
    
    # set the format of edge
    E(g)$label <- E(g)$text
    E(g)$label.cex <- 0.8
    E(g)$arrow.size <- 0.3
    
    # plot the tree
    plot(g, layout = layout_as_tree)
}


plotHeatMap <- function(accu) {
    library(ggplot2)
    dimnames(accu) <- list(1:nrow(accu),1:ncol(accu))
    accu.m <- melt(accu)
    
    for (i in 1:ncol(accu)) {
        accu.m[(i - 1)*nrow(accu) + (1:nrow(accu)),] <- accu.m[(i-1)*nrow(accu) + (nrow(accu):1),]
    }
    
    names(accu.m) <- c('nmin','minleaf','err')
    print(accu.m)
    
    p <- ggplot(accu.m, aes(minleaf,nmin)) + 
        geom_tile(aes(fill = err), colour = "white") +
        scale_fill_gradient(low = "green", high = "red")
    return(p)
}

plotP3D <- function(accu) {
    library(plot3D)
    persp3D(z = accu, theta = 50,phi = 25, expand = 0.75, 
            xlab = 'nmin',ylab = 'minleaf',zlab = 'err',axes = TRUE)
}

