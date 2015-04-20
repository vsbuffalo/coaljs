
function range(from, to, by) {
    var by = typeof(by) === 'undefined' ? 1 : by;
    var out = [];
    for (var i = from; i <= to; i += by) out.push(i);
    return out;
}

function countLeaves(tree) {
    var count = 0;
    function leafCounter(subtree) {
	if (typeof subtree.children === 'undefined') {
	    count++;
	    return;
	}
	subtree.children.forEach(function(l) {
	    leafCounter(l);
	});
    }
    leafCounter(tree);
    return count;
}

function depth(depth) {
    var tmp = function(x) { return x.depth == depth; };
    return tmp;
}

function GeneologyTree(tree, options) {
    if (typeof options === 'undefined')
	options = {height: 440, width: 440}; // TODO
    // functions for recursing through a tree (nested object) and grabbing all
    // edges and nodes.
    // Coal time parameters
    var t_time = 0;
    function total_time(tree) {
	if (tree.children) {
	    t_time+= tree.time;
	    tree.children.forEach(total_time);
	}
    };
    total_time(tree); // get total time by recursing tree
    
    function nodes() {
	// Get all nodes in a nice array. This has side effects: it will append
	// slots.
	var _nodes = [], depth = 0, center = 0.5;
	var get_nodes = function(node, slots, is_left) {
	    node.is_left = is_left;
	    node.slots = slots.slice(0);
	    node.depth = depth;
	    node.y = typeof node.time === 'undefined' ? 0 : node.time / t_time;
	    _nodes.push(node);
	    if (node.children) {
		depth++;
		if (node.children.length != 2)
		    throw new Error("length children != 2");
		// get the nodes for left and right children
		var left = node.children[0]; 
		var right = node.children[1]; 
		center = center/2;
		left.x = node.x - center;
		right.x = node.x + center;
		get_nodes(left, slots.splice(0, countLeaves(left)), true);
		get_nodes(right, slots.splice(0, countLeaves(right)), false);
	    }
	};
	tree.x = 0.5;
	get_nodes(tree, range(0, countLeaves(tree)-1), null);
	return _nodes;
    };

    function links() {
	// get all edges as a link (parent |---| child) for all nodes, return in
	// arary
	var _edges = [];
	var get_edges = function(node, depth) {
	    var depth = typeof depth  == 'undefined' ? 0 : depth + 1;
	    if (typeof node.children != 'undefined') {
		node.children.forEach(function(b) {
		    var edge = {
			parent: node,
			child: b,
			depth: depth
		    };
		    _edges.push(edge);
		    get_edges(b, depth);
		});
	    }
	};
	get_edges(tree);
	return _edges;
    };
    return {total_time: t_time, tree: tree, nodes: nodes, links: links};
}





