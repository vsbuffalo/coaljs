
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

function CoalescentGeneology(k) {
    var geno = [], samples = range(0, k-1).map(function(o) {return {name: o};});
    var coal, coal_time = 0, coal_times =
	    range(2, k).map(function(i) { return Dist.exp(i*(i-1)/2); });
    // draw two lineages two coalesce, and pop up a coalescent time for them
    while (coal_times.length) {
	coal_time += coal_times.pop()
	coal = {time: coal_time,
		children: sampleArrayWithoutReplacement(samples, 2)};
	samples.push(coal);
    }
    return samples[0];
}

/// Test
var margin = {top: 40, right: 10, bottom: 20, left: 10},
    width = 400 - margin.right - margin.left,
    height = 500 - margin.top - margin.bottom;

var svg = d3.select("body").append("svg")
	.attr("width", width + margin.right + margin.left)
	.attr("height", height + margin.top + margin.bottom);


// Test
var coal = CoalescentGeneology(6),
    nodes = GeneologyTree(coal).nodes();
var tree = d3.layout.tree();

function xScale(x) {
    return width*x;
}

function yScale(y) {
    return height*y;
}


svg.selectAll("circles")
    .data(GeneologyTree(coal).nodes())
    .enter().append("circle")
    .attr("id", function(o) { return "ct" + o.time; })
    .attr("cx", function(o) { return xScale(o.x); })
    .attr("cy", function(o) { return height - yScale(o.y); }) 
    .attr("r", "4")
    .style("fill", function(o) {
	if (o.children)
	    return "steelblue";
	return "black";
    });
	

svg.selectAll("lines")
    .data(GeneologyTree(coal).links())
    .enter().append("line")
	.attr({
	    x1: function(n) { return xScale(n.parent.x); },
	    x2: function(n) { return xScale(n.child.x); },
	    y1: function(n) { return height - yScale(n.parent.y); },
	    y2: function(n) { return height - yScale(n.child.y); }
	})
	.style("fill", "none")
	.style("stroke", "gray")
	.style("stroke-opacity", 0.4);











