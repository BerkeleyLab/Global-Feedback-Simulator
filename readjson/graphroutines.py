import pydot
def filter_dotify(fil,graphin,fname="test.png"):
    graph = pydot.Dot(graph_type='digraph',rankdir="LR")
    innode = pydot.Node("Input")
    graph.add_node(innode)
    las = innode
    for o in xrange(fil.order):
        if fil.A_modes[o]>1:
            summer = pydot.Node("Sigma{0}".format(o),label="&Sigma;",shape="circle")
            graph.add_node(summer)

        for m in xrange(fil.A_modes[o]):
            nod = pydot.Node(str(fil.A_poles[fil.A_coeff_start[o]+m]))
            graph.add_node(nod)
            #these += [nod]
            if fil.A_modes[o]>1:
                graph.add_edge(pydot.Edge(nod,summer))
            edg = pydot.Edge(las,nod)
            graph.add_edge(edg)
        if fil.A_modes[o]>1:
            las=summer
        else:
            las=nod

    ounode = pydot.Node("Output")
    graph.add_node(ounode)
    graph.add_edge(pydot.Edge(las,ounode))
    graph.write_png(fname)
