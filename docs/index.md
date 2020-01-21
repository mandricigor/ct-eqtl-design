

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
		<td colspan="2">
			<p><div id="buttonPreset1"></div> We aim to have a cell-type specific eQTL study (ct-eQTL) for B cells with available budget $35,000. As the optimal coverage is 10,000 reads per cell (<a href="https://www.biorxiv.org/content/10.1101/766972v1">Mandric et al, 2019</a>), what is the sample size and the number of cells we should have per individual? And what will be the effective sample size?</p>
		</td>
	</tr>
	<tr valign="bottom">
		<td>
                        Cell type
			<div id="dropCelltype" style="width:250px;"></div>
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
		<td>
			Budget ($1,000)
			<div id="sliderBudget" style="width:250px;"></div><input type="text" id="inpBudget" value="15" />
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
<div style="margin:0 auto;"><img id="catImg" src="" height="0"></img></div>
</div>
<div id="plot" style="width:689px;height:350px;position:absolute;"></div>
</div>

<table id="inputTable">
<tr>
<td colspan="2">
<div id="description" style="width:689px;position:absolute;">
<p>
		This calculator provides guidance in selecting optimal experimental designs for cell-type specific eQTL analysis (ct-eQTL). It is based on the results of the paper "Optimal design of single-cell RNA sequencing experiments for cell-type specific eQTL analysis" <a href="https://www.biorxiv.org/content/10.1101/766972v1">(Mandric et al, 2019)</a>. As it is shown in their paper, the optimal coverage at which one should conduct cell-type specific eQTL studies is 10,000 (±2,500) reads per cell. The assumptions underlying this calculator are the following:
<ul>
  			<li>Genetic multiplexing is used (At most 16 samples per pool, at most 24,000 cells per pool) </li>
  			<li>The cost of library prep is $2,000 per reaction</li>
  			<li>The cost of sequencing is $5 per million reads</li>
</ul> 
</p>
<p>
	This calculator was developed by Igor Mandric (UCLA) and Christoph Hafemeister (NYU).
</p>
</div>


<script>
// min, max, step, default
//var sampleSizeRange = [10, 1000, 1, 60];
var budgetRange = [5, 100, 1, 15];

var slopes = {"CD14+ Monocytes": 1.556, "CD4 T cells": 1.549, "B cells": 1.197, "CD8 T cells": 1.196, "NK cells": 1.239, "Megakaryocytes": 1.273, "FCGR3A+ cells": 1.304, "Dendritic cells": 1.180};
var intercepts = {"CD14+ Monocytes": 4.846, "CD4 T cells": 5.074, "B cells": 4.645, "CD8 T cells": 1.623, "NK cells": 2.763, "Megakaryocytes": 1.103, "FCGR3A+ cells": 3.614, "Dendritic cells": 2.875};


$("#cat").hide(0);

$("#inpBudget").val(budgetRange[3]);
//$("#inpSampleSize").val(sampleSizeRange[3]);
$("#celltype :selected").val("NK cells");

function loadPreset(values) {
	$("#inpBudget").val(values[0]);
	$("#inpSampleSize").val(values[1]);
        $('#celltype').val(values[2]);
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
            if (info.length == batch_ind_info.length) {
                // group by batch cell count
                info_singlet_reads_dict = {};
                info_multiplet_reads_dict = {};
                info_singlets_pic_dict = {};
                info_nonident_multi_pic_dict = {};
                //for entry in info:
                for (i = 0; i < info.length; i ++) {
                    entry = info[i];
                    info_singlets_pic_dict[entry[2]] = getget(info_singlets_pic_dict, entry[2]) + entry[0];
                    info_nonident_multi_pic_dict[entry[3]] = getget(info_nonident_multi_pic_dict, entry[3]) + entry[0];
                    info_singlet_reads_dict[entry[4]] = getget(info_singlet_reads_dict, entry[4]) + entry[0];
                    info_multiplet_reads_dict[entry[5]] = getget(info_multiplet_reads_dict,entry[5]) + entry[0];
                }
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
                reads_pp.push([cn, singlets_pic, nonident_multi_pic, singlets_reads, multiplets_reads]);
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


function optimal_designs(budget) {
    var lowCell = 500;
    var highCell = 2750;
    var lowInd = 5;
    var highInd = 1000;
    var uu = exp_design_fixed_lane_capacity(budget, lowCell, highCell, lowInd, highInd);
    good_ind = {};
    for (const u in uu) {
        var uarr = uu[u];
        var i;
        for (i = 0; i < uarr.length; i ++) {
            if ((uarr[i][3] > 7500) && (uarr[i][3] < 12500)) {
                //console.log(u, uarr[i][0]);
                if (u in good_ind) {
                    if (uarr[i][0] > good_ind[i]) {
                        good_ind[i] = uarr[i][0];
                    }
                }
                else {
                    good_ind[u] = uarr[i][0];
                }
            }
        }
    }
    console.log(good_ind);
    return good_ind;
}


function updateResults() {
    syncInput();
    myslope = slopes[$('#celltype :selected').text()];
    myintercept = intercepts[$('#celltype :selected').text()];
    mybudget = $('#inpBudget').val();

    var optimal = optimal_designs(parseInt(mybudget) * 1000);

    var ess = mybudget * myslope + myintercept;

    $("#results").empty();	
    if (Object.keys(optimal).length > 0) {
        $("#results").append('<p>We recommend the following experimental designs:</p>');
        $("#results").append("<ul>");
        for (const ii in optimal) {
            $("#results").append("<li>" + ii + " individuals and " + optimal[ii] + " cells per individual</li>");
        }
        $("#results").append("</ul>");
    }
    else {
        $("#results").append('<p>Sorry, there was a problem</p>');
    }
    //$("#results").append('<p>Effective Sample Size (ESS): '+Plotly.d3.format(",.r")(ess.toFixed(0))+'</p>');

    var x = [];
    var y = [];
    var x2 = [mybudget, mybudget];
    var y2 = [0 * myslope + myintercept, mybudget * myslope + myintercept];
	for (var i = 5; i <= 100; i += 1) {
		x.push(i);
		y.push(i * myslope + myintercept);
	}

	// for the plot
	var c1 = 'rgb(27,158,119)';
	var c2 = 'rgb(217,95,2)';
	var layout = {
		margin: {l: 60, r: 60, b: 60, t: 10, pad: 0},
		showlegend: false,
		xaxis: {title: 'Budget ($1,000)', zeroline: false, fixedrange: true},
		yaxis: {title: 'Effective Sample Size',
				titlefont: {color: c1},
    			tickfont: {color: c1},
    			hoverformat: '.0f',
    			zeroline: false,
    			fixedrange: true},
                annotations: [{x:mybudget, y:ess, xref:"x", yref:"y", text:"ESS:" + ess.toFixed(0), align:"center", bordercolor: '#c7c7c7',borderwidth: 2,bgcolor: '#ff7f0e',opacity: 0.8, arrowcolor: '#636363',arrowwidth: 2,arrowsize: 1,arrowhead: 2, font:{size:20, color:"#ffffff"}}]
	}
	var trace1 = {x: x, y: y, line:{color: c1}, name: ''};
	var trace2 = {x: x2, y: y2, line:{color: c2, width:3, dash: "dot", mode:"lines"}, hoverinfo:"skip"};
        var trace3 = {x:[mybudget], y: [mybudget * myslope + myintercept], mode: "markers", marker:{color: 'rgb(142, 124, 195)',size:20}, name: "ESS", hoverinfo:"skip"};
	var plotElem = document.getElementById('plot');
	Plotly.purge(plotElem);
	Plotly.plot(plotElem, [trace1, trace2, trace3], layout, {displayModeBar: false});
	if (catMode > 0 && !$("#cat").is(":visible")) {
		showCat();
		setTimeout(hideCat, 1000);
	}
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

//$("#sliderSampleSize").slider({
//	value: parseInt($("#inpSampleSize").val()),
//	min: sampleSizeRange[0],
//	max: sampleSizeRange[1],
//	step: sampleSizeRange[2],
//	slide: function(event, ui) {$("#inpSampleSize").val(ui.value); updateResults();},
//	stop: function(event, ui) {updateResults();},
//});


$("#buttonPreset1").button({label: 'Load preset'}).click(function() {loadPreset([35, 60, "B cells"]);});


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

