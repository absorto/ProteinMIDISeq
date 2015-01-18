from nodebox.graphics import *
from nodebox.graphics.physics import Node, Edge, Graph

# Create a graph with randomly connected nodes.
# Nodes and edges can be styled with fill, stroke, strokewidth parameters.
# Each node displays its id as a text label, stored as a Text object in Node.text.
# To hide the node label, set the text parameter to None.
abo = list("ATGGCCGAGGTGTTGCGGACGCTGGCCGGAAAACCAAAATGCCACGCACTTCGACCTATGATCCTTTTCCTAATAATGCTTGTCTTGGTCTTGTTTGGTTACGGGGTCCTAAGCCCCAGAAGTCTAATGCCAGGAAGCCTGGAACGGGGGTTC")

g = Graph()
# bases as nodes
index = 1
for b in abo:
    g.add_node(id=index, 
               radius = 6,
               stroke = color(0), 
               text = color(.3,.6,.9))
    index += 1
# Random edges.
for n in range(0,len(abo)-1):
    node1 = g.nodes[n]
    node1.text = Text(abo[n], font="Droid Serif", fontsize=5, fontweight=BOLD) 
    node2 = g.nodes[n+1]
    g.add_edge(node1, node2, 
        length = 1.0, 
        weight = 1.0, 
        stroke = color(.7,.1,.1))

#g.prune(depth=0)          # Remove orphaned nodes with no connections.
g.distance         = 5   # Overall spacing between nodes.
g.layout.force     = 0.01 # Strength of the attractive & repulsive force.
g.layout.repulsion = 15   # Repulsion radius.

dragged = None
def draw(canvas):
    
    canvas.clear()
    background(1)
    translate(250, 250)
    
    # With directed=True, edges have an arrowhead indicating the direction of the connection.
    # With weighted=True, Node.centrality is indicated by a shadow under high-traffic nodes.
    # With weighted=0.0-1.0, indicates nodes whose centrality > the given threshold.
    # This requires some extra calculations.
    g.draw(weighted=False, directed=True)
    g.update(iterations=2)
    
    # Make it interactive!
    # When the mouse is pressed, remember on which node.
    # Drag this node around when the mouse is moved.
    dx = canvas.mouse.x - 250 # Undo translate().
    dy = canvas.mouse.y - 250
    global dragged
    if canvas.mouse.pressed and not dragged:
        dragged = g.node_at(dx, dy)
    if not canvas.mouse.pressed:
        dragged = None
    if dragged:
        dragged.x = dx
        dragged.y = dy
        
canvas.size = 1280, 800
canvas.run(draw)
