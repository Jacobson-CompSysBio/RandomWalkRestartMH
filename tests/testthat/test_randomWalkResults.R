library(RandomWalkRestartMH)

############################################################
context("Random Walk Computation")
############################################################

## Multiplex
m1 <- igraph::make_graph(c(1,2,1,3,2,3), directed = FALSE)
m2 <- igraph::make_graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
multiObject_1 <- create.multiplex(list(m1=m1,m2=m2))
TranMatrix <- compute.transition.matrix.homogeneous(multiObject_1)
Multiplex1_Seeds <- c(1)

RWR_MultiResults <- 
    Random.Walk.Restart.Multiplex(TranMatrix, multiObject_1, Multiplex1_Seeds)

RWRM_ExpectedResults <- data.frame(NodeNames = c("3","2","4"),
                                    Score = c(0.040624865, 
                                    0.0161811209, 0.0135287952),
                                    stringsAsFactors = FALSE)
print(RWR_MultiResults)
## Multiplex-Heterogeneous
h1 <- igraph::make_graph(c("A","C","B","E","E","D","E","C"), directed = FALSE)
multiObject_2 <- create.multiplex(list(h1=h1))
bipartite_relations <- data.frame(m=c(1,3),h=c("A","E"))

multiHetObject <- 
    create.multiplexHet(multiObject_1, multiObject_2, bipartite_relations)
MultiHetTranMatrix <- compute.transition.matrix.heterogeneous(multiHetObject)

Multiplex2_Seeds <- c("E")
 
RWR_MultiHetResults <- Random.Walk.Restart.MultiplexHet(MultiHetTranMatrix,
    multiHetObject,Multiplex1_Seeds,Multiplex2_Seeds)

print(RWR_MultiHetResults)
RWRMH_Multi1ExpectedResults <- data.frame(NodeNames = c("3","2","4"),
                                          Score = c(0.03948193089443327,
                                                    0.004671625367690084,
                                                    0.003384466823594482), 
                                         stringsAsFactors = FALSE)

RWRMH_Multi2ExpectedResults <- data.frame(NodeNames = c("A","C","B","D"),
                                          Score = c(0.06316305892641791,
                                                    0.028345414073428166,
                                                    0.018870955342772815,
                                                    0.018870955342772815),
                                          stringsAsFactors = FALSE)


rownames(RWR_MultiHetResults$RWRMH_Multiplex1) <- NULL
rownames(RWRMH_Multi1ExpectedResults) <- NULL

rownames(RWR_MultiHetResults$RWRMH_Multiplex2) <- NULL
rownames(RWRMH_Multi2ExpectedResults) <- NULL

test_that("Random Walk Results are correct", {
    expect_equal(RWR_MultiResults$RWRM_Results, RWRM_ExpectedResults, 
                tolerance = 0.00001)

    # TODO: These tests test heterogeneous multiplex layers. As we are currently in discussion on
    # how to move forward with the construction of these layers, we are not expending energy on
    # updating these test cases.  Please reach out to @kpsmithjr or @lanematthewj for questions
    # expect_equal(RWR_MultiHetResults$RWRMH_Multiplex1,
    #             RWRMH_Multi1ExpectedResults,
    #             tolerance = 0.00001)
    # expect_equal(RWR_MultiHetResults$RWRMH_Multiplex2, 
    #             RWRMH_Multi2ExpectedResults,
    #             tolerance = 0.00001)
})
