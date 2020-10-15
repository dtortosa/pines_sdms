function (tree, x, res = 100, fsize = NULL, ftype = NULL, lwd = 4, 
    legend = NULL, lims = NULL, outline = TRUE, sig = 3, type = "phylogram", 
    direction = "rightwards", plot = TRUE, ...) 
{

tree = tree_prunned
x = seed_mass_vector
res = 100
fsize = NULL
ftype = NULL
lwd = 4
legend = NULL
lims = NULL
outline = TRUE
sig = 3
type = "phylogram"
direction = "rightwards"
plot=FALSE
method = "fastAnc"


    if (!inherits(tree, "phylo")){ 
        stop("tree should be an object of class \\"phylo\\".")
    }    
    if (hasArg(mar)){
        mar <- list(...)$mar
    } else {
        mar <- rep(0.3, 4)
    }
    if (hasArg(offset)){
        offset <- list(...)$offset       
    } else {
        offset <- NULL
    }
    if (hasArg(method)){
        method <- list(...)$method
    } else {
        method <- "fastAnc"
    }
    if (hasArg(hold)){
        hold <- list(...)$hold
    } else {
        hold <- TRUE
    }
    if (hasArg(leg.txt)){ 
        leg.txt <- list(...)$leg.txt
    } else {
        leg.txt <- "trait value"
    }

    #altura máxima de un nodo (nodos más alejado de la raiz o finales de rama -> actualidad)
    h <- max(nodeHeights(tree))

    #set the gradient color in basis on the maximun height (nodo most separated) and the resoution of the gradient (res)
    steps <- 0:res/res * max(h)

    #calcular la altura sobre la raiz de cada nodo (como de lejos están de la raíz)
    #se calcula para los dos nodos de cada rama, así todas las ramas con 85.73146 son ramas terminales, porque 85.73146 es la altura máxima
    H <- nodeHeights(tree)

    #reconstruction of ancestral state of nodes 
    if (method == "fastAnc"){ #fastAnc
        a <- fastAnc(tree, x)
    } else{
        if (method == "anc.ML") { #anc.ML
            fit <- anc.ML(tree, x)
            a <- fit$ace
        }    
        if (!is.null(fit$missing.x)){ 
            x <- c(x, fit$missing.x)
        }     
    } else { 
        if (method == "user"){ #el usuario provee los estados ancestrales en anc.states
            if (hasArg(anc.states)){ #está anc.states usado como argumento en la función?
                anc.states <- list(...)$anc.states
            } else {
                cat("No ancestral states have been provided. Using states estimated with fastAnc.\\n\\n")
                a <- fastAnc(tree, x)
            }
        
            #si los estados que aportaas son menos que el número total de nodos entonces hayq eu recosnturir los que faltan. 
            if (length(anc.states) < tree$Nnode) {
                nodes <- as.numeric(names(anc.states))
                tt <- tree
                for (i in 1:length(nodes)) {
                    M <- matchNodes(tt, tree, method = "distances")
                    ii <- M[which(M[, 2] == nodes[i]), 1]
                    tt <- bind.tip(tt, nodes[i], edge.length = 0, 
                      where = ii)
                }
                a <- fastAnc(tt, c(x, anc.states))
                M <- matchNodes(tree, tt, method = "distances")
                a <- a[as.character(M[, 2])]
                names(a) <- M[, 1]
            }
            else {
                if (is.null(names(anc.states))) 
                    names(anc.states) <- 1:tree$Nnode + Ntip(tree)
                a <- anc.states[as.character(1:tree$Nnode + Ntip(tree))]
            }
        }
    }

    y <- c(a, x[tree$tip.label])
    names(y)[1:Ntip(tree) + tree$Nnode] <- 1:Ntip(tree)
    A <- matrix(y[as.character(tree$edge)], nrow(tree$edge), ncol(tree$edge))
    cols <- rainbow(1001, start = 0, end = 0.7)
    names(cols) <- 0:1000
    if (is.null(lims)){ 
        lims <- c(min(y), max(y))
    }    
    trans <- 0:1000/1000 * (lims[2] - lims[1]) + lims[1]
    names(trans) <- 0:1000
    tree$maps <- vector(mode = "list", length = nrow(tree$edge))

    #for each rama
    for (i in 1:nrow(tree$edge)) {
        XX <- cbind(
            c(H[i, 1], #altura del primer nodo de la rama (el más cercano a la raiz ó más antiguo)
            steps[intersect(which(steps > H[i, 1]), which(steps < H[i, 2]))]), 
            c(steps[intersect(which(steps > H[i, 1]), which(steps < H[i, 2]))], 
            H[i, 2])) - H[i, 1]
        YY <- rowMeans(XX)
        if (!all(YY == 0)) {
            b <- vector()
            for (j in 1:length(YY)){
                b[j] <- (A[i, 1]/YY[j] + 
                A[i, 2]/(max(XX) - YY[j]))/(1/YY[j] + 1/(max(XX) - 
                YY[j]))
            }    
        } else {
            b <- A[i, 1]
        }
        d <- sapply(b, getState, trans = trans)
        tree$maps[[i]] <- XX[, 2] - XX[, 1]
        names(tree$maps[[i]]) <- d
    }
    tree$mapped.edge <- makeMappedEdge(tree$edge, tree$maps)
    tree$mapped.edge <- tree$mapped.edge[, order(as.numeric(colnames(tree$mapped.edge)))]
    class(tree) <- c("simmap", setdiff(class(tree), "simmap"))
    xx <- list(tree = tree, cols = cols, lims = lims)
    class(xx) <- "contMap"
    if (plot) 
        plot.contMap(xx, fsize = fsize, ftype = ftype, lwd = lwd, 
            legend = legend, outline = outline, sig = sig, type = type, 
            mar = mar, direction = direction, offset = offset, 
            hold = hold, leg.txt = leg.txt)
    invisible(xx)
}