var COLORS = ['#8b0000','#91000a','#960212','#9d0519','#a2081f','#a70d25',
	    '#ac122c','#b11732','#b61b37','#ba213e','#bf2543','#c32b48',
	    '#c7304d','#cb3453','#cf3a58','#d23f5d','#d64563','#da4967',
	    '#dd4e6c','#e05471','#e35976','#e65e7a','#e9647f','#eb6984',
	    '#ee6e88','#f0748c','#f27990','#f47e94','#f68599','#f88a9c',
	    '#fa8fa1','#fb95a4','#fd9ba8','#fea0ac','#ffa6b0','#ffacb3',
	    '#ffb1b7','#ffb7bb','#ffbdbe','#ffc4c2','#ffc9c5','#ffd0c9',
	    '#ffd4cb','#ffdbce','#ffe1d2','#ffe5d4','#ffebd7','#fff2da',
	    '#fff6dc','#fffcdf','#fbffdd','#f3fdd7','#ecfcd2','#e4facc',
	    '#def8c8','#d8f6c4','#d1f3c0','#ccf2bd','#c6efba','#c0edb6',
	    '#bbebb3','#b7e8b1','#b0e6ad','#ace4ab','#a7e1a9','#a2dea7',
	    '#9edba4','#98d9a2','#95d6a1','#90d39f','#8bd19d','#88ce9b',
	    '#83cb9a','#7fc998','#7ac696','#77c495','#73c094','#6fbd92',
	    '#6bbb91','#67b890','#64b58f','#5fb28e','#5baf8d','#57ac8c',
	    '#53aa8b','#50a78a','#4ca589','#48a288','#449e87','#3f9c87',
	    '#3c9986','#389685','#339384','#2e9184','#288e83','#248b82',
	    '#1e8882','#168581','#0c8281','#008080'];


function range(from, to, by) {
    var by = typeof(by) === 'undefined' ? 1 : by;
    var out = [];
    for (var i = from; i <= to; i += by) out.push(i);
    return out;
}

function sampleArrayWithoutReplacement(x, n) {
    n = typeof x === 'undefined' ? 1 : n;
    var sample = [];
    if (n > x.length) throw new Error("n must be <= length of x");
    while (n) {
	sample.push(x.splice(Dist.dunif(0, x.length-1), 1)[0]);
	n--;
    }
   return sample;
}

function CoalescentGenealogy(k, theta) {
    var lineages = range(0, k-1).map(function(o) {return {name: o, mutations: []};});
    var total_time = 0, mutations = [], event_time,
	is_mutation, i, nsites = 0;
    var cols = COLORS.slice(0); // duplicate for this sim
    while (lineages.length > 1) {
	i = lineages.length;
	event_time = Dist.exp(i*(i-1)/2 + i*theta/2);
	total_time += event_time;
	is_mutation = Dist.bern(theta/(i + theta - 1));
	if (is_mutation) {
	    // draw random lineage to put mutation on
	    var lineage = lineages[Dist.dunif(0, lineages.length-1)];
	    var mutation = {lineage: lineage,
			    total_time: total_time,
			    time: event_time,
			    color: sampleArrayWithoutReplacement(cols, 1),
			    id: "mut-"+nsites,
			    site: nsites};
	    lineage.mutations.push(mutation);
	    mutations.push(mutation);
	    nsites++;
	} else {
	    // is coalescent event
	    var coalescing_lineages = sampleArrayWithoutReplacement(lineages, 2);
	    var coalesce = {lineages: coalescing_lineages,
			    time: event_time,
			    mutations: [],
			    total_time: total_time};
	    coalescing_lineages.forEach(function(o) { return o.parent = coalesce;});
	    lineages.push(coalesce);
	}
    };
    return {lineages: lineages[0],
	    mutations: mutations};
}

function GenealogyTree(coal, options) {
    var tree = coal.lineages;
    var muts = coal.mutations;
    if (typeof options === 'undefined')
	options = {stem:0.1}; // TODO
    // functions for recursing through a tree (nested object) and grabbing all
    // edges and nodes.
    // Coal time parameters
    var t_time = tree.total_time;
    var max_depth = 0, depth = 0, left_depth, right_depth;
    var leaves = [];
    function maxDepth(node) {
	if (node.lineages) {
	    node.lineages[0].depth = node.depth + 1;
	    node.lineages[1].depth = node.depth + 1;
	    left_depth = maxDepth(node.lineages[0]);
	    right_depth = maxDepth(node.lineages[1]);
	    if (left_depth > right_depth)
		return left_depth + 1;
	    else
		return right_depth + 1;
	} else {
	    leaves.push(node);
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
	    node.y = typeof node.time==='undefined' ? 0 : node.total_time*0.3;
	    if (typeof node.y === 'undefined') debugger;
	    _nodes.push(node);
	    if (node.lineages) {
		if (node.lineages.length != 2)
		    throw new Error("length lineages != 2");
		// get the nodes for left and right lineages
		var left = node.lineages[0]; 
		var right = node.lineages[1]; 
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

    function mutations() {
	return muts.map(function(o) {
	    o.x = o.lineage.x;
	    return o;
	});
    };

    function links() {
	// get all edges as a link (parent |---| child) for all nodes, return in
	// arary
	var _edges = [];
	var get_edges = function(node, depth) {
	    var depth = typeof depth  == 'undefined' ? 0 : depth + 1;
	    if (typeof node.lineages != 'undefined') {
		node.lineages.forEach(function(b) {
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

    function propagate_leaf_mutations() {
	function add_mutations(node, leaf) {
	    if (node.parent) {
		Array.prototype.push.apply(leaf.mutations, node.parent.mutations);
		add_mutations(node.parent, leaf);
	    } 
	};
	leaves.forEach(function(o) { add_mutations(o, o); });
    };
    propagate_leaf_mutations();

    
    return {max_depth: max_depth, total_time: t_time, tree: tree,
	    mutations: mutations, leaves: leaves,
	    nodes: nodes, links: links};
}


/// Test
var margin = {top: 40, right: 10, bottom: 20, left: 10},
    width = 800 - margin.right - margin.left,
    height = 600 - margin.top - margin.bottom;

var tree_height = 0.75*height;

function pi(tree) {
    // TODO: worth testing this a bit more
    var leaves = tree.leaves,
	muts = tree.mutations(), pi = 0;
    for (var i = 0; i < leaves.length; i++) {
	for (var j = 0; j < i; j++) {
	    var diff = leaves[j].mutations.filter(function(x) { return leaves[i].mutations.indexOf(x) < 0; }).length + 
		    leaves[i].mutations.filter(function(x) { return leaves[j].mutations.indexOf(x) < 0; }).length;
	    pi += (1/(leaves.length-1))*(1/leaves.length)*diff;
	}
    }
    return 2*pi;
};

function newSim() {
    var svg = d3.select("body").append("svg")
	    .attr("id", "coalsim")
	    .attr("width", width + margin.right + margin.left)
	    .attr("height", height + margin.top + margin.bottom);


    var coal = CoalescentGenealogy(9, 10);
    var tree = GenealogyTree(coal);
    var nodes = tree.nodes();
    var links = tree.links();
    var muts = tree.mutations();
    var leaves = tree.leaves;

    function xScale(x) {
	return width*x;
    }

    function yScale(y) {
	return tree_height*y;
    }

    svg.selectAll("circles")
	.data(nodes)
	.enter().append("circle")
	.attr("id", function(o) { return "ct" + o.coal_time; })
	.attr("cx", function(o) { return xScale(o.x); })
	.attr("cy", function(o) { return tree_height - yScale(o.y); }) 
	.attr("r", 1)
	.style("fill", "gray");

    svg.selectAll("polyline")
	.data(links)
	.enter().append("polyline")
	.attr("points", function(n) {
	    var points = [[xScale(n.parent.x), tree_height-yScale(n.parent.y)],
			  [xScale(n.child.x), tree_height-yScale(n.parent.y)],
			  [xScale(n.child.x), tree_height-yScale(n.child.y)]];
	    return points.map(function(x) { return x.join(",") }).join(" ");
	})
	.style("fill", "none")
	.style("stroke-width", 2)
	.style("stroke", "#CCC");

    svg.selectAll("text")
	.data(nodes)
	.enter().append("text")
	.filter(function(x) { return typeof x.lineages === 'undefined'; })
	.attr("y", function(n) { return 16 + tree_height-yScale(n.y); })
	.attr("x", function(n) { return -4 + xScale(n.x); })
	.text(function(n) { return n.name; })
	.style("fill", "gray");

    svg.selectAll("circles")
	.data(muts)
	.enter().append("circle")
	.attr("class", "mut")
	.attr("id", function(o) { return o.id; })
	.attr("cx", function(o) { return xScale(o.x); })
	.attr("cy", function(o) { return tree_height - yScale(o.total_time*0.3); }) 
	.attr("r", 3)
	.style("fill", function(o) { return o.color;});


    var nmuts = muts.length;

    svg.selectAll("sequences")
	.data(leaves)
	.enter().append("polyline")
	.attr("points", function(n) {
	    var y = height - 10*(n.name);
	    n.y = y;
	    // while we're here, process mutations
	    n.mutations = n.mutations.map(function(m) {
     		m.name = n.name; // the leaf name
		if (typeof m.y === "undefined")
		    m.y = [];
		m.y.push(y);
		m.x =  width/10 + (width*0.8/muts.length)*m.site;
     		return m;
     	    });
	    var points = [[width/10, y], [width - width/10, y]];
	    return points.map(function(x) { return x.join(",") }).join(" ");
	})
	.style("stroke-width", 1)
	.style("stroke", "#CCC");

    svg.selectAll("sequence-labels")
	.data(leaves)
	.enter().append("text")
	.attr("x", function(n) { return width/10 - 8; })
	.attr("y", function(n) { return n.y + 4; })
	.text(function(n) { return n.name; })
	.style({"font-size": 10,
		"fill": "#CCC"});




    // a big of data munging to make new mutations for each y value
    // this should probably be done in a d3 command.
    var all_leaf_muts = [];
    leaves.forEach(function(n) {
	n.mutations.forEach(function(m) {
	    m.y.forEach(function(y) {
		all_leaf_muts.push({y: y, x: m.x, id: m.id, color:m.color});
	    });
	});
    });

    svg.selectAll("leaf-seqs")
	.data(all_leaf_muts)
	.enter().append("circle")
	.attr("class", "mut")
	.attr("id", function(n) { return n.id; })
	.attr("cy", function(n) { return n.y; })
	.attr("cx", function(n) { return n.x; })
	.attr("r", 3)
	.style("fill", function(o) { return o.color;});

    var mutalles = svg.selectAll(".mut");

    mutalles.on("mouseover", function(x) {
	d3.selectAll("#"+x.id)
	    .attr("r", 5);
    });

    mutalles.on("mouseout", function(x) {
	d3.selectAll("#"+x.id)
	    .attr("r", 3);
    });

    d3.selectAll("#stats")
	.text(function() {
	    var a = d3.sum(range(1, leaves.length-1).map(function(x) { return 1/x; }));
	    var str = "summary statistics: S = " + muts.length + ", θ = " + Math.round(100*muts.length/a)/100 + ", π = " + Math.round(100*pi(tree))/100;
	    return str;
	});
    return tree;
}

newSim();


d3.select("#refresh").on("click", function() {
    d3.select("#coalsim").remove();
    // new simulation
    newSim();
});

// var pis = [];
// for (var i = 0; i < 100; i++) {
//     var tree= newSim();
//     pis.push(pi(tree));
// }

// console.log(d3.mean(pis));
