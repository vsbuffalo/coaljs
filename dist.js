
var m = new MersenneTwister();

var Dist = (function(unif) {
    // discrete uniform
    function dunif(min, max) {
	return Math.floor(Math.random() * (max - min + 1)) + min;
    };

    // Poisson rv
    function pois(lambda) {
	var n = 0,
	    limit = Math.exp(-lambda),
	    x = m.random();
	while (x > limit) {
	    n++;
	    x *= m.random();
	}
	return n;
    };

    // Exponential rv
    function exp(lambda) {
	return -Math.log(m.random())/lambda;
    };

    function bern(p) {
	return Number(m.random() < p);
    };

    // binomial rv
    function binom(n, p) {
	var x = 0;
	for (var i = 0; i < n; i++) x += bern(p);
	return x;
    };
    return {dunif: dunif,
	    pois: pois,
	    exp: exp,
	    bern: bern,
	    binom: binom};
})();
