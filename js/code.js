function $S(selector) { return document.querySelector(selector); }
function $Sall(selector) { return document.querySelectorAll(selector); }

const {
    ChiSquared,
    Multinomial,
    rng: { MersenneTwister },
    rng: { normal: { Inversion } },
    R: { map, sum, div, mult, multiplex, numberPrecision, seq: _seq }
} = libR;

const mt = new MersenneTwister();
const { dmultinom, rmultinom } = Multinomial(mt);

const { qchisq } = ChiSquared();

const pow = multiplex(Math.pow);
const add = multiplex((x,e)=>e+x);
const sumrow = multiplex(sum);
const addrow = function(arr, x) { return map(arr)((e)=>add(e,x)); }
const divrow = function(arr, x) { return map(arr)((e)=>div(e,x)); }
// multiplex doesn't handle arr / arr operations nicely
const powrow = multiplex(pow);

var verbose, confidence, Ntests, precision, N, p, ptest, dim, chisq_sig;

function parsePage() {
    verbose = $S('#verbose').checked;
    confidence = Number($S('#confidence').value);
    Ntests = Number($S('#Ntests').value);
    precision = Number($S('#precision').value);
    N = Number($S('#Nsamples').value);
    [_,p] = parseArray($S('#p').value);
    [_,ptest] = parseArray($S('#ptest').value);
    [goodToGo,dim] = initCheck(p, ptest);
    if (goodToGo) {
	chisq_sig = qchisq(confidence, dim-1);
    }
    return goodToGo;
}

function initCheck(p, ptest) {
    let dim = p.length;
    if (ptest.length != dim) {
	addStr("ERROR: p and ptest are not the same length");
	return [false, NaN];
    }
    return [true, dim];
}

function accuracyFromSample(N,solo) {
    if (!N) { N = window.N; }
    let counts = rmultinom(Ntests,N,ptest);
    // chisq_arr = N*np.sum((counts/(1.0*N) - p)**2/p, axis=1)
    let tmp1a = divrow(counts, N);
    let tmp1b = mult(p,-1);
    let tmp2 = addrow(tmp1a, tmp1b);
    let tmp3 = powrow(tmp2, 2);
    let tmp4 = divrow(tmp3, p);
    let tmp5 = sumrow(tmp4);
    let tmp6 = mult(tmp5, N);
    let chisq_arr = tmp6;
    let sigCount = chisq_arr.filter((e)=>e>chisq_sig).length;
    let acc = sigCount/Ntests
    if (verbose || solo) {
	addStr(`significance of ${N} = ${acc}${solo ? ` w/ ${100*confidence}% confidence.` : ''}`);
    }
    return acc;
}

function sampleForConfidence() {
    let N = 1;
    let accuracy = 0.0;
    let delta = 1e10;
    let initIncrease = true;
    let isIncrease = false;
    let swapsAtOne = 0;
    while ((Math.abs(accuracy - confidence) > precision || accuracy < confidence)) {
        if (initIncrease) {
            if (accuracy > confidence) {
                delta = N/4.0;
                N = N-delta;
                initIncrease = false;
                if (verbose) {
		    addStr(`inital increase finished, d=${delta}`);
		}
            } else { N = N*2; }
        } else {
            if ((isIncrease && accuracy > confidence) ||
                (!isIncrease && accuracy < confidence)) {
                delta = Math.max(1, delta/2.0);
		if (delta<=1) {
		    if (swapsAtOne > 2) {
			addStr(`search got stuck around N=${N} and accuracy=${accuracy} for confidence in [${confidence}, ${confidence+precision}], try running again with widened precision`);
			break;
		    }
		    swapsAtOne += 1;
		}
                if (verbose) {
		    addStr(`passed confidence, d=${delta}`);
		}
                isIncrease = !isIncrease;
	    }
            N = N + (isIncrease ? 1.0 :  -1.0) * delta;
	}
	accuracy = accuracyFromSample(N);
    }
    if (swapsAtOne <= 2) {
	addStr(`Number of Tests: ${Ntests}\nTo differentiate between p [${p}] and ptest [${ptest}] with ${100*confidence}% confidence, use at least\nSample Size: ${N} (accuracy = ${accuracy})`);
    }
    return N;
}

function addStr(strtmp) {
    $S('#output').innerText += strtmp + '\n';
}

function parseFraction(val) {
    if (!val.includes('/')) return NaN
    return val.split('/').reduce((p,c)=>p/c);
} //  https://stackoverflow.com/questions/7142657/

// attempts to coerce val into a number
// tries  convert to float and convert to fraction
// Number is stricter than parseFloat
function parseVal(val) {
    return parseFraction(val) || Number(val);
}

function isValid(val, numval) {
    return val!="" && !Number.isNaN(numval)
}

function parseArray(arrstr) {
    let regex = /\[(?:[^,]+,)*(?:[^,])+\]/
    // starts with [ ends with ]
    // contains chars separated by , last digit is not ,
    if (!regex.test(arrstr)) {
	return [false, "ERROR: Incorrect array format, must be [p1,p2,...,pN"];
    }
    let strarr = arrstr.substr(1,arrstr.length-2).split(',');
    let numarr = [];
    for (let strnum of strarr) {
	let numnum = parseVal(strnum);
	if (!isValid(strnum, numnum)) {
	    return [false, `ERROR: ${strnum} is not a number (in ${arrstr})`];
	}
	numarr.push(numnum)
    }
    let totalp = sum(numarr);
    if (totalp <= 0) {
        return [false, `ERROR: total probability should be greater than 0, not ${totalp}`];
    }
    if (totalp != 1) {
        numarr <- numarr / totalp
        addStr(`Total probability != 1. Stardardizing probability vector to [${numarr}]`);
    }
    return [true, numarr]
}

function clearOutput() {
    $S('#output').innerText = '';
}

function updateNumberInput(e) {
    let id = e.srcElement.id;
    let val = e.srcElement.value;
    let num = Number(val);
    if (! isValid(val, num)) {
	e.srcElement.value = 1;
	e.srcElement.dispatchEvent(new Event('change'));
	addStr(`ERROR: ${id} '${val}' is not a number`);
	return false;
    } else {
	let min = e.srcElement.min!='' ? Number(e.srcElement.min) : -Infinity;
	let max = e.srcElement.max!='' ? Number(e.srcElement.max) : Infinity;
	if (num < min) {
	    e.srcElement.value = min;
	    e.srcElement.dispatchEvent(new Event('change'))
	    addStr(`ERROR: '${val}' < minimum ${min} for ${id}`);
	    return false;
	} else if (num > max) {
	    e.srcElement.value = max;
	    e.srcElement.dispatchEvent(new Event('change'))
	    addStr(`ERROR: '${val}' > maximum ${max} for ${id}`);
	    return false;
	}
    }
    return true;
}

function updateConfidenceText() {
    let tmpc = Number($S('#confidence').value);
    let tmpp = Number($S('#precision').value);
    $S('#sampleForN').innerText = `find N with confidence in [${tmpc}, ${Math.min(1.0,tmpc+tmpp)}]`;
}

window.addEventListener('load', function() {
    updateConfidenceText();
    for (let node of $Sall('input[type=number]')) {
	node.addEventListener('change', (e)=>{
	    if (updateNumberInput(e)) {
		let id = e.srcElement.id
		if (id=='confidence'||id=='precision') {
		    updateConfidenceText();
		}
	    }
	});
    }
    for (let node of $Sall('tr[name=input2d]>td>input')) {
	node.addEventListener('change', (e)=>{
	    let val = e.srcElement.value;
	    let pval = parseVal(val);
	    let elID = e.srcElement.id;
	    let tgtID = elID=="test_p" ? "#ptest" : "#p"
	    if (! isValid(val, pval)) {
		addStr(`ERROR: ${elID} '${val}' is not a number`);
		e.srcElement.value = 1;
		e.srcElement.dispatchEvent(new Event('change'));
	    } else if (pval<0 || pval>1) {
		addStr(`ERROR: ${elID} '${pval}' must be within [0,1]`);
		e.srcElement.value = pval<0 ? 0 : 1
		e.srcElement.dispatchEvent(new Event('change'));
	    } else {
		$S(tgtID).value = `[${pval},${1-pval}]`;
		$S(tgtID).dispatchEvent(new Event('change'));
	    }
	});
    }
    $S('#outputClear').addEventListener('click', clearOutput);
    for (let node of $Sall('input[name=pdim]')) {
	node.addEventListener('change', (e)=>{
	    if (e.srcElement.value=='2D') {
		// switch TO 2D
		for (let node of $Sall('tr[name=input2d]')) {
		    node.classList.remove('hide')
		}
		for (let node of $Sall('input[name=pvals]')) {
		    node.disabled = true
		}
	    } else {
		// switch FROM 2D
		for (let node of $Sall('tr[name=input2d]')) {
		    node.classList.add('hide')
		}
		for (let node of $Sall('input[name=pvals]')) {
		    node.disabled = false
		}
	    }
	});
    }
    for (let node of $Sall('tr[name=arrnd]>td>input')) {
	node.addEventListener('change', (e)=>{
	    let ID = e.srcElement.id;
	    let [isValid, output] = parseArray(e.srcElement.value);
	    if (!isValid) {
		addStr(output);
		let tgtID = ID=='p' ? '#null_p' : '#test_p';
		$S(tgtID).dispatchEvent(new Event('change'));
	    }
	});
    }
    for (let node of $Sall('tr[name=execButton]>td>button')) {
	node.addEventListener('click', (e)=>{
	    if (parsePage()) {
		clearOutput();
		let id = e.srcElement.id;
		if (verbose || 'accFromN') {
		    addStr(`chisq significance threshold = ${chisq_sig} with df=${dim-1} and ${100*confidence}% confidence`);
		}
		if (id=='accFromN') {
		    accuracyFromSample(N,true);
		} else if (id=='sampleForN') {
		    sampleForConfidence();
		} else {
		    console.error('this should not happen');
		}
	    }
	});
    }
});
