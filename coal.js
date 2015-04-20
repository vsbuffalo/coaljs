
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
    var coal, coal_times =
	    range(2, k).map(function(i) { return Dist.exp(i*(i-1)/2); });
    // draw two lineages two coalesce, and pop up a coalescent time for them
    while (coal_times.length) {
	coal = {time: coal_times.pop(),
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
var coal = GeneologyTree(CoalescentGeneology(10)).nodes();
var tree = d3.layout.tree();
svg.selectAll("circles")
    .data(GeneologyTree(coal).nodes())
    .enter().append("circle")
    .attr("cx", function(o) { return 20 + 100*o.x; })
    .attr("cy", function(o) { return 20 + (1-o.time)*400; }) 
    .attr("r", "4")
    .style("fill", function(o) {
	if (o.children)
	    return "steelblue";
	return "black";
    });
	



// svg.selectAll("circles")
//     .data(GeneologyTree(coal).nodes())
//     .enter().append("circle")
//     .attr("cx", function(o) { return 20 + 100*o.x; })
//     .attr("cy", function(o) { return 20 + (1-o.time)*400; }) 
//     .attr("r", "4")
//     .style("fill", function(o) {
// 	if (o.children)
// 	    return "steelblue";
// 	return "black";
//     });
	
