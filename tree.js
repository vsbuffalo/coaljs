
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

function GenealogyTree(tree, options) {
    if (typeof options === 'undefined')
	options = {stem:0.1}; // TODO
    // functions for recursing through a tree (nested object) and grabbing all
    // edges and nodes.
    // Coal time parameters
    var t_time = tree.total_time;
    var max_depth = 0, depth = 0, left_depth, right_depth;
    function maxDepth(node) {
	if (node.children) {
	    node.children[0].depth = node.depth + 1;
	    node.children[1].depth = node.depth + 1;
	    left_depth = maxDepth(node.children[0]);
	    right_depth = maxDepth(node.children[1]);
    if (left_depth > right_depth)
		return left_depth + 1;
	    else
		return right_depth + 1;
	} else {
	    return 0;
	}
    };
    tree.depth = 1;
    max_depth = maxDepth(tree); // get total time by recursing tree

    function nodes() {
	// Get all nodes in a nice array. This has side effects: it will append
	// slots.
	var _nodes = [], depth = 0, center = 0.5/max_depth,
	    ymax = t_time/(1-options.stem);
	var get_nodes = function(node, slots, is_left) {
	    node.is_left = is_left;
	    node.slots = slots.slice(0);
	    node.y = typeof node.coal_time==='undefined' ? 0 : node.total_time*0.3;
	    _nodes.push(node);
	    if (node.children) {
		if (node.children.length != 2)
		    throw new Error("length children != 2");
		// get the nodes for left and right children
		var left = node.children[0]; 
		var right = node.children[1]; 
		var left_slots = slots.splice(0, countLeaves(left));
		var right_slots = slots.splice(0, countLeaves(right));
		left.x = node.x - center/node.depth;
		right.x = node.x + center/node.depth;
		get_nodes(left, left_slots, true);
		get_nodes(right, right_slots, false);
	    }
	    //console.log(center);
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
    return {max_depth: max_depth, total_time: t_time, tree: tree, nodes: nodes, links: links};
}





