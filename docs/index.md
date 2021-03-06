

  [Calculator](index.md) |
  [About](about.md)



<head>
    <link rel="stylesheet" href="https://code.jquery.com/ui/1.12.1/themes/base/jquery-ui.css">
	<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
	<script src="https://code.jquery.com/jquery-3.1.1.min.js"></script>
	<script src="https://code.jquery.com/ui/1.12.1/jquery-ui.min.js"></script>
	<script src="https://underscorejs.org/underscore-min.js"></script>
	
</head>

<table id="inputTable">
	<tr>
		<td colspan="3">
			<p><div id="buttonPreset1"></div> We aim to have a cell-type specific eQTL study (ct-eQTL) for B cells with available budget $35,000. As the optimal coverage at an average sequencing cost of 5$ per 1 million reads is 10,000 reads per cell (<a href="https://www.biorxiv.org/content/10.1101/766972v1">Mandric et al, 2019</a>), what is the sample size and the number of cells we should have per individual? And what will be the effective sample size?</p>
		</td>
	</tr>
	<tr valign="bottom">
		<td style="width:33%;">
                        Cell type
			<div id="dropCelltype"></div>
    				<select name="celltype" id="celltype">
      					<option selected="selected">CD4 T cells</option>
      					<option>CD14+ Monocytes</option>
					<option>B cells</option>
					<option>CD8 T cells</option>
      					<option>NK cells</option>
      					<option>FCGR3A+ cells</option>
					<option>Megakaryocytes</option>
      					<option>Dendritic cells</option>
    				</select>
		</td>
		<td style="width:34%;">
			Budget ($1,000)
			<div id="sliderBudget"></div><input type="text" id="inpBudget" value="55" />
		</td>
		<td style="width:33%;">
			Sequencing cost (per 1M reads)
			<div id="dropSeqCost"></div>
				<select name="seqcost" id="seqcost">
					<option>1 $</option>
					<option>2 $</option>
					<option>3 $</option>
					<option>4 $</option>
					<option selected="selected">5 $</option>
				</select>
		</td>
	</tr>
</table>



<table width="100%">
<tr>
<td valign="top" colspan="2">
<div id="results" style="width:689px;"></div>
</td>
</tr></table>

<div style="width:689px;height:350px;">
<div id="cat" style="width:689px;height:350px;position:absolute;display:flex;flex-direction:column;">
<div style="margin:0 auto;">Calculating results</div>
</div>
<div id="plot" style="width:689px;height:350px;position:absolute;"></div>
</div>

<table id="inputTable">
<tr>
<td colspan="2">
<div id="description" style="width:689px;position:absolute;">
<p>
	<br>
	<h3>Main assumptions</h3>
	<p style="text-align:left">
		This calculator provides guidance in selecting optimal experimental designs for cell-type specific eQTL analysis (ct-eQTL). It is based on the results of the paper "Optimal design of single-cell RNA sequencing experiments for cell-type specific eQTL analysis" <a href="https://www.biorxiv.org/content/10.1101/766972v1">(Mandric et al, 2019)</a>. As it is shown in their paper, the optimal coverage at which one should conduct cell-type specific eQTL studies is 10,000 (±2,500) reads per cell. The assumptions underlying this calculator are the following:
	</p>
<ul style="text-align:left">
  			<li>Genetic multiplexing is used</li>
                            <ul>
                                <li>At most 16 samples per pool</li>
                                <li>At most 24,000 cells per pool</li>
                            </ul>
  			<li>The cost of library prep is $2,000 per reaction</li>
  			<li>The cost of sequencing is $5 per million reads</li>
</ul> 
</p>
<p>
	<h3>How to interpret the output of the calculator?</h3>
    <ul style="text-align:left">
        <li>The area of best experimental designs is marked with a blue dashed line</li>
        <li>The Optimal Effective Sample Size is the same for all the experimental designs inside the marked area</li>
        <li>The blue dots show some particular experimental designs. Hover the mouse upon them to see the number of pools and number of reads per cell</li>
        <li>The number of cells stands for <i>desired</i> number of cells</li>
    </ul>
</p>
<p>
	This calculator was developed by <a href="mailto:imandric@ucla.edu" target="_top">Igor Mandric (UCLA)</a> and Christoph Hafemeister (NYGC).
</p>
</div>


<script>
// min, max, step, default
//var sampleSizeRange = [10, 1000, 1, 60];
var budgetRange = [10, 100, 1, 55];

var slopes = {"CD14+ Monocytes": 1.556, "CD4 T cells": 1.549, "B cells": 1.197, "CD8 T cells": 1.196, "NK cells": 1.239, "Megakaryocytes": 1.273, "FCGR3A+ cells": 1.304, "Dendritic cells": 1.180};
var intercepts = {"CD14+ Monocytes": 4.846, "CD4 T cells": 5.074, "B cells": 4.645, "CD8 T cells": 1.623, "NK cells": 2.763, "Megakaryocytes": 1.103, "FCGR3A+ cells": 3.614, "Dendritic cells": 2.875};
var magnitude = {"1 $": 1.75, "2 $": 1.5, "3 $": 1.2, "4 $": 1.1, "5 $": 1.0}
var lows = {"1 $": 15000, "2 $": 12500, "3 $": 10500, "4 $": 9500, "5 $": 7500}
var highs = {"1 $": 20000, "2 $": 17500, "3 $": 15500, "4 $": 14500, "5 $": 12500}

$("#cat").hide(0);

$("#inpBudget").val(budgetRange[3]);
//$("#inpSampleSize").val(sampleSizeRange[3]);
$("#celltype :selected").val("NK cells");
$("#seqcost :selected").val("5 $");

function loadPreset(values) {
	$("#inpBudget").val(values[0]);
	$("#inpSampleSize").val(values[1]);
        $('#celltype').val(values[2]);
	$('#seqcost').val(values[3]);
	updateResults();
}

function checkInput(range, txtInp, sliderInp) {
	var value = parseFloat($(txtInp).val());
	if (isNaN(value)) value = range[3];
	if (value < range[0]) value = range[0];
	if (value > range[1]) value = range[1];
	$(txtInp).val(value);
	$(sliderInp).slider("option", "value", value);
}


function syncInput() {
	checkInput(budgetRange, "#inpBudget", "#sliderBudget");
	//checkInput(sampleSizeRange, "#inpSampleSize", "#sliderSampleSize");
        $("#celltype").selectmenu("refresh");
	$("#seqcost").selectmenu("refresh");
}


LIBRARY_PREP_COST = 2000
ILLUMINA_PER_MILLION = 5
MULTIFACTOR = 1.82
R = 0.5714286
M = 4.5997701e-6


function p(multi) {
    return multi / (M * (1 - multi));
}


function q(nr, multi){
    return - (nr * multi) / (R * M * (1 - multi));
}

function numCellsLoaded(cells, multi){
    return -0.5 * p(multi) - Math.sqrt(0.25 * p(multi) * p(multi) - q(cells, multi));
}

function multiplet_rate(ncl){
    return M * ncl;
}


function numCellsRecovered(cells, multi) {
    return R * numCellsLoaded(cells, multi);
}


function singlet_rate(ncl){
    return 1 - multiplet_rate(ncl);
}


function num_singlet(cells, multi) {
    return parseInt(singlet_rate(numCellsLoaded(cells, multi)) * numCellsRecovered(cells, multi));
}

function num_ident_multiplet(cells, multi){
    numMultiplet = numCellsRecovered(cells, multi) - num_singlet(cells, multi);
    return numMultiplet * (multi - 1) / multi;
}


function num_multiplet(cells, multi){
    return numCellsRecovered(cells, multi) - num_singlet(cells, multi);
}

function num_nonident_multiplet(cells, multi){
    return parseInt(num_multiplet(cells, multi) - num_ident_multiplet(cells, multi));
}


function readsX(cells, reads_pc, multi){
    nsing = num_singlet(cells, multi);
    nmult = num_multiplet(cells, multi);
    nidentmulti = num_nonident_multiplet(cells, multi);
    return parseInt(cells * reads_pc) / ((nsing / (nsing + MULTIFACTOR * nmult)) + (nidentmulti / (1/MULTIFACTOR * nsing + nmult)));
}


function singletAvgReadsX(cells, reads_pc, multi){
    rx = readsX(cells, reads_pc, multi);
    nsing = num_singlet(cells, multi);
    nmulti = num_multiplet(cells, multi);
    return parseInt(rx / (nsing + MULTIFACTOR * nmulti));
}

function multiAvgReadsX(cells, reads_pc, multi) {
    rx = readsX(cells, reads_pc, multi);
    nsing = num_singlet(cells, multi);
    nmulti = num_multiplet(cells, multi);
    return parseInt(rx / ((1/MULTIFACTOR) * nsing + nmulti));
}



function dichotomy(cells, money, multi, eps=0.00001){
    mini = 1
    maxi = 1000000
    var i;
    for (i=0; i < 20; i ++){
        midi = parseInt(0.5 * (mini + maxi));
        midi_reads = readsX(cells, midi, multi);
        money2 = midi_reads * ILLUMINA_PER_MILLION / 1000000;
        if (Math.abs((money2 - money) * 1.0 / money) < eps){
            return midi;
        }
        else if (money2 > money) {
            maxi = midi;
        }
        else if (money2 <= money) {
            mini = midi;
        }
    }
    return midi
}


function exp_design(budget, lo_cell, hi_cell, lo_p, hi_p, diff_cell=250, multi=8) {
    // ASSUMPTION 1: number of persons is divisible by multi(=16)
    // ASSUMPTION 2: number of cells is divisible by diff_cell(=250)
    design = {};
    pers = hi_p;
    while (pers >= lo_p) {
        seq_budget = budget - (pers / multi) * LIBRARY_PREP_COST;
        if (seq_budget < 0) {
            design[pers] = new Array();
            break;
        }
        // find budget per sequencing batch
        seq_batch_budget = seq_budget / (pers / multi);
        reads_pp = new Array();
        cn = hi_cell;
        while (cn >= lo_cell) {
            cells_batch = cn * multi;
            rpp = dichotomy(cells_batch, seq_batch_budget, multi);
            if (rpp > 0) {
                singlets_ = num_singlet(cells_batch, multi);
                nonident_multi_ = num_nonident_multiplet(cells_batch, multi);
                singlets_reads = singletAvgReadsX(cells_batch, rpp, multi);
                multiplets_reads = multiAvgReadsX(cells_batch, rpp, multi);
                singlets_pic = parseInt(singlets_ / multi);
                nonident_multi_pic = parseInt(nonident_multi_ / multi);
                reads_pp.push([cn, singlets_pic, nonident_multi_pic, singlets_reads, multiplets_reads]);
            cn -= diff_cell;
            }
        }
        if (reads_pp) {
            design[pers] = reads_pp;
        }
        pers -= multi;
    }
    return design
}


//var uu = exp_design(35000, 500, 2750, 40, 120);


function sum(arr){
  return arr.reduce(function(a,b){
    return a + b
  }, 0);
}

function getget(myObj, el) {
    if (el in myObj) {
        return myObj[el];
    }
    else {
        return 0;
    }
}

function exp_design_fixed_lane_capacity(budget, lo_cell, hi_cell, lo_p, hi_p, diff_cell=250, diff_person=8, capacity=24000, max_multi=16) {
    // ASSUMPTION 1: number of cells per lane is maximized
    // ASSUMPTION 2: maximum number of individuals multiplexed is 16
    // Put greedily cells into lanes taking care to not exceed the maximum lane capacity
    // and not to exceed number of multiplexed persons
    design = {};
    pers = hi_p;
    while (pers >= lo_p) {
        cn = hi_cell;
        reads_pp = new Array();
        while (cn >= lo_cell) {
            number_ind_per_lane = parseInt(capacity * 1.0 / cn);
            number_ind_per_lane = Math.min(number_ind_per_lane, max_multi);
            nr_batches = parseInt(pers * 1.0 / number_ind_per_lane);
            if (pers % number_ind_per_lane > 0) {
                nr_batches += 1;
            }
            number_ind_per_lane_approx = pers / nr_batches;
            total_seq_budget = budget - LIBRARY_PREP_COST * nr_batches;
            if (total_seq_budget <= 0) {
                break;
            }
            seq_budget_per_person = total_seq_budget / pers;
            batch_ind_info = new Array();
            var i;
            for (i = 0; i < nr_batches; i ++) {
                batch_ind_info.push(number_ind_per_lane_approx);
            }
            extras = pers - sum(batch_ind_info);
            //cyc = cycle(range(nr_batches)) # STOPPED HERE
            cyc = 0;
            while (extras > 0) {
                inc_batch = cyc % nr_batches;
                batch_ind_info[inc_batch] += 1;
                extras -= 1;
                cyc += 1;
            }
            batch_money_info = new Array();
            for (i = 0; i < batch_ind_info.length; i ++) {
                batch_money_info.push(seq_budget_per_person * batch_ind_info[i]);
            }
            info = new Array();
            for (i = 0; i < batch_ind_info.length; i ++) {
                v = batch_money_info[i];
                w = batch_ind_info[i];
                u = cn * w;
                rpp = dichotomy(u, v, w);
                if (rpp > 0) {
                    singlets_ = num_singlet(u, w);
                    nonident_multi_ = num_nonident_multiplet(u, w);
                    singlets_reads = singletAvgReadsX(u, rpp, w);
                    multiplets_reads = multiAvgReadsX(u, rpp, w);
                    singlets_pic = parseInt(singlets_ / w);
                    nonident_multi_pic = parseInt(nonident_multi_ / w);
                    info.push([w, cn, singlets_pic, nonident_multi_pic, singlets_reads, multiplets_reads]);
                }
            }
            //console.log(batch_ind_info); console.log("VAFLI");
            //console.log(info); console.log("SUKASUKA");
            if (info.length == batch_ind_info.length) {
                // group by batch cell count
                info_singlet_reads_dict = {};
                info_multiplet_reads_dict = {};
                info_singlets_pic_dict = {};
                info_nonident_multi_pic_dict = {};
                //for entry in info:
                for (i = 0; i < info.length; i ++) {
                    entry = info[i];
                    //console.log(entry);
                    //console.log("ENTRY");
                    info_singlets_pic_dict[entry[2]] = getget(info_singlets_pic_dict, entry[2]) + entry[0];
                    info_nonident_multi_pic_dict[entry[3]] = getget(info_nonident_multi_pic_dict, entry[3]) + entry[0];
                    info_singlet_reads_dict[entry[4]] = getget(info_singlet_reads_dict, entry[4]) + entry[0];
                    info_multiplet_reads_dict[entry[5]] = getget(info_multiplet_reads_dict,entry[5]) + entry[0];
                }
                //console.log(info_singlets_pic_dict);
                //console.log(info_singlet_reads_dict);
                singlets_pic = 0;
                singlets_pic_sum = 0;
                for (const u in info_singlets_pic_dict) {
                    v = info_singlets_pic_dict[u];
                    singlets_pic += u * v;
                    singlets_pic_sum += v;
                }
                singlets_pic /= parseFloat(singlets_pic_sum);
                singlets_pic = parseInt(singlets_pic);
                nonident_multi_pic = 0;
                nonident_multi_pic_sum = 0;
                for (const u in info_nonident_multi_pic_dict) {
                    v = info_nonident_multi_pic_dict[u];
                    nonident_multi_pic += u * v;
                    nonident_multi_pic_sum += v;
                }
                nonident_multi_pic /= parseFloat(nonident_multi_pic_sum);
                nonident_multi_pic = parseInt(nonident_multi_pic);
                singlets_reads = 0;
                singlets_reads_sum = 0;
                for (const u in info_singlet_reads_dict) {
                    v = info_singlet_reads_dict[u];
                    singlets_reads += u * v;
                    singlets_reads_sum += v;
                }
                singlets_reads /= parseFloat(singlets_reads_sum);
                singlets_reads = parseInt(singlets_reads);
                multiplets_reads = 0;
                multiplets_reads_sum = 0;
                for (const u in info_multiplet_reads_dict) {
                    v = info_multiplet_reads_dict;
                    multiplets_reads += u * v;
                    multiplets_reads_sum += v;
                }
                multiplets_reads /= parseFloat(multiplets_reads_sum);
                multiplets_reads = parseInt(multiplets_reads);
                reads_pp.push([cn, singlets_pic, nonident_multi_pic, singlets_reads, multiplets_reads, info.length]);
            }
            cn -= diff_cell;
            if (reads_pp) {
                design[pers] = reads_pp;
            }
        }
        pers -= diff_person;
    }
    return design;
}


function optimal_designs(budget, low_cov, high_cov) {
    var lowCell = 500;
    var highCell = 2750;
    var lowInd = 10;
    var highInd = 1000;
    var uu = exp_design_fixed_lane_capacity(budget, lowCell, highCell, lowInd, highInd);
    //console.log(uu);
    //console.log("vasea");
    good_ind = {};
    best_designs = [];
    for (const u in uu) {
        good_ind[u] = [];
        var uarr = uu[u];
        var i;
        for (i = 0; i < uarr.length; i ++) {
            if ((uarr[i][3] > low_cov) && (uarr[i][3] < high_cov)) {
                console.log(u, uarr[i]);
                //if (u in good_ind) {
                //    if (uarr[i][0] > good_ind[i]) {
                //        good_ind[i] = [uarr[i][0], uarr[i][5], Math.round( (uarr[i][3] / 1000) * 10 ) / 10];
                //    }
                //}
                //else {
                    best_designs.push([parseInt(u), uarr[i][0], uarr[i][5], Math.round( (uarr[i][3] / 1000) * 10 ) / 10]);
                //}
            }
        }
    }
    //console.log(good_ind);
    return best_designs;
    //return good_ind;
}

function formatNumber(num) {
  return num.toString().replace(/(\d)(?=(\d{3})+(?!\d))/g, '$1,')
}


var convexhull = new function() {
	
	// Returns a new array of points representing the convex hull of
	// the given set of points. The convex hull excludes collinear points.
	// This algorithm runs in O(n log n) time.
	this.makeHull = function(points) {
		var newPoints = points.slice();
		newPoints.sort(this.POINT_COMPARATOR);
		return this.makeHullPresorted(newPoints);
	};
	
	
	// Returns the convex hull, assuming that each points[i] <= points[i + 1]. Runs in O(n) time.
	this.makeHullPresorted = function(points) {
		if (points.length <= 1)
			return points.slice();
		
		// Andrew's monotone chain algorithm. Positive y coordinates correspond to "up"
		// as per the mathematical convention, instead of "down" as per the computer
		// graphics convention. This doesn't affect the correctness of the result.
		
		var upperHull = [];
		for (var i = 0; i < points.length; i++) {
			var p = points[i];
			while (upperHull.length >= 2) {
				var q = upperHull[upperHull.length - 1];
				var r = upperHull[upperHull.length - 2];
				if ((q.x - r.x) * (p.y - r.y) >= (q.y - r.y) * (p.x - r.x))
					upperHull.pop();
				else
					break;
			}
			upperHull.push(p);
		}
		upperHull.pop();
		
		var lowerHull = [];
		for (var i = points.length - 1; i >= 0; i--) {
			var p = points[i];
			while (lowerHull.length >= 2) {
				var q = lowerHull[lowerHull.length - 1];
				var r = lowerHull[lowerHull.length - 2];
				if ((q.x - r.x) * (p.y - r.y) >= (q.y - r.y) * (p.x - r.x))
					lowerHull.pop();
				else
					break;
			}
			lowerHull.push(p);
		}
		lowerHull.pop();
		
		if (upperHull.length == 1 && lowerHull.length == 1 && upperHull[0].x == lowerHull[0].x && upperHull[0].y == lowerHull[0].y)
			return upperHull;
		else
			return upperHull.concat(lowerHull);
	};
	
	
	this.POINT_COMPARATOR = function(a, b) {
		if (a.x < b.x)
			return -1;
		else if (a.x > b.x)
			return +1;
		else if (a.y < b.y)
			return -1;
		else if (a.y > b.y)
			return +1;
		else
			return 0;
	};
	
};




function updateResults() {
    syncInput();
    myslope = slopes[$('#celltype :selected').text()];
    myintercept = intercepts[$('#celltype :selected').text()];
    mybudget = $('#inpBudget').val();
    myseqcost = $('#seqcost').val();

    var optimal = optimal_designs(parseInt(mybudget) * 1000, lows[myseqcost], highs[myseqcost]);
    console.log(optimal);

    console.log("AAAA");
    console.log(magnitude[myseqcost]);

    var ess = (mybudget * myslope + myintercept) * magnitude[myseqcost];

    $("#results").empty();	
    //if (1 == 0) {
    //if (Object.keys(optimal).length > 0) {
    //    $("#results").append('<p>We recommend the following experimental designs:</p>');
    //    $("#results").append("<ul>");
    //    for (const ii in optimal) {
    //        $("#results").append("<li>" + ii + " individuals and " + optimal[ii] + " cells per individual</li>");
    //    }
    //    $("#results").append("</ul>");
    //}
    //else {
    //    $("#results").append('<p>Sorry, there was a problem</p>');
    //}
    //$("#results").append('<p>Effective Sample Size (ESS): '+Plotly.d3.format(",.r")(ess.toFixed(0))+'</p>');

    var points = [];
    var inds = [];
    var cells = [];
    var pools = [];
    var coverages = []; // reads per cell
    var hovers = [];
    for (i = 0; i < optimal.length; i ++) {
        inds.push(optimal[i][0]);
        cells.push(optimal[i][1]);
        points.push({x: optimal[i][0], y: optimal[i][1]});
        pools.push(optimal[i][2]);
        coverages.push(optimal[i][3] * 1000);
        hovers.push(sprintf("Reads per cell: %s<br>Number of pools: %s", formatNumber(optimal[i][3] * 1000), optimal[i][2]));
    }

    console.log(inds);
    console.log(cells);
    console.log(pools);
    console.log(coverages);
    console.log(points);
    console.log("GHGHGHGHGHG");

    var hull = convexhull.makeHull(points);
    var hullx = [];
    var hully = [];
    for (i = 0; i < hull.length; i ++) {
        hullx.push(hull[i].x);
        hully.push(hull[i].y);
    }
    hullx.push(hull[0].x);
    hully.push(hull[0].y);

    var trace11 = {
        x: inds,
        y: cells,
        text: hovers,
        mode: 'markers',
        type: 'scatter',
        marker: {size: 15},
        hovertemplate:
            "<b>%{text}</b><br><br>" +
            "%{yaxis.title.text}: %{y:,.0f}<br>" +
            "%{xaxis.title.text}: %{x:,.0f}<br>" +
            "<extra></extra>"
    };

    trace22 = {
        x: hullx,
        y: hully,
        mode: 'lines',
        line: {dash: "dashdot", color: "blue", width: 1}
    }

    var maxind = Math.max.apply(null, inds) + 10;
    var minind = Math.min.apply(null, inds) - 10;
    var maxcell = Math.max.apply(null, cells) + 300;
    var mincell = Math.min.apply(null, cells) - 300;

    data = [trace11, trace22];
    var layout = {
        showlegend: false,
        hovermode: "closest",
        xaxis: {
        tickwidth: 4,
        ticklen: 8,
        range: [minind, maxind],
        tickfont: {
            size: 20,
            color: "green"
        },
        title: "Individuals",
        titlefont: {
            size: 25
        }
      },
      yaxis: {
        tickwidth: 4,
        ticklen: 8,
        range: [mincell, maxcell],
        tickfont: {
            size: 13,
            color: "green"
        },
        title: "Cells",
        titlefont: {
            size: 25
        }
      },
      title: {
          text: sprintf('Optimal Effective Sample Size: %s', Math.round(ess, 0)),
          font: {
              size: 24
          }
      }
    };
    //Plotly.newPlot('results', data, layout);
    var plotElem = document.getElementById('plot');
    Plotly.purge(plotElem);
    Plotly.plot(plotElem, data, layout, {displayModeBar: false});
}




$("#inputTable").find("td").css("padding", "12px");

$("#celltype").selectmenu({
	value: $("#celltype").val(),
        change: function(event, ui) {updateResults();},
});

$("#sliderBudget").slider({
	value: parseInt($("#inpBudget").val()),
	min: budgetRange[0],
	max: budgetRange[1],
	step: budgetRange[2],
	slide: function(event, ui) {$("#inpBudget").val(ui.value); updateResults();},
	stop: function(event, ui) {updateResults();},
});

$("#seqcost").selectmenu({
	value: $("#seqcost").val(),
        change: function(event, ui) {updateResults();},
});

//$("#sliderSampleSize").slider({
//	value: parseInt($("#inpSampleSize").val()),
//	min: sampleSizeRange[0],
//	max: sampleSizeRange[1],
//	step: sampleSizeRange[2],
//	slide: function(event, ui) {$("#inpSampleSize").val(ui.value); updateResults();},
//	stop: function(event, ui) {updateResults();},
//});


$("#buttonPreset1").button({label: 'Load preset'}).click(function() {loadPreset([35, 60, "B cells", "5 $"]);});


$("#moreDetails").button({label: 'More detail'});
$("#moreDetails").click(function() {
	$(this).text(function(i, text){
    	return text === "More detail" ? "Less detail" : "More detail";
    });
    updateResults();
});

$(".ui-button").css('padding', 2);

var lazyUpdate = _.debounce(updateResults, 500);
$("input").keyup(lazyUpdate);

updateResults();

</script>


