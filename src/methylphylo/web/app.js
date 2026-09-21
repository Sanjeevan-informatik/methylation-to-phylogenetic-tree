'use strict';
const $ = id => document.getElementById(id);
let recipes = {}, engineAvailable = false, inputText = '', inputType = '', zipUrl, svgUrl;
const examples = {
  dna: '>A\nACGTACGTACGTACGTACGT\n>B\nACGTACGTACGTACGTACGA\n>C\nATGTATGTATGTACGTACGT\n>D\nATGTATGTATGTACGTACGA\n',
  beta: 'locus\tA\tB\tC\tD\ncg1\t0.1\t0.12\t0.8\t0.82\ncg2\t0.2\t0.22\t0.7\t0.75\ncg3\t0.15\t0.13\t0.9\t0.88\ncg4\t0.8\t0.85\t0.1\t0.12\ncg5\t0.5\t0.52\t0.48\t0.5\n',
  binary: 'locus\tA\tB\tC\tD\ncg1\t0\t0\t1\t1\ncg2\t0\t0\t1\t1\ncg3\t1\t1\t0\t0\ncg4\t0\t1\t0\t1\ncg5\t0\t0\t0\t0\n',
  'binary-fasta': '>A\n00100\n>B\n00110\n>C\n11000\n>D\n11010\n',
  vcf: '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD\n1\t100\t.\tA\tG\t60\tPASS\t.\tGT\t0/0\t0/0\t1/1\t1/1\n1\t200\t.\tC\tT\t60\tPASS\t.\tGT\t0/0\t0/0\t1/1\t1/1\n1\t300\t.\tG\tA\t60\tPASS\t.\tGT\t0/0\t1/1\t0/0\t1/1\n'
};
examples.counts = 'locus\tsample\tmethylated\ttotal\n' + ['cg1','cg2','cg3','cg4'].flatMap((l,i)=>['A','B','C','D'].map((s,j)=>`${l}\t${s}\t${j<2?i+1:9-i}\t10\n`)).join('');
function updateRecipes() {
  const type = $('data-type').value;
  const kind = ['dna','vcf'].includes(type) ? 'dna' : ['beta','counts'].includes(type) ? 'beta' : 'binary';
  $('recipe').replaceChildren();
  for (const [key,r] of Object.entries(recipes)) {
    if (r.kind !== kind) continue;
    const option = new Option(r.label + (r.method === 'ml' && !engineAvailable ? ' · IQ-TREE required' : ''),key);
    option.disabled = r.method === 'ml' && !engineAvailable;
    $('recipe').add(option);
  }
  $('engine').textContent = engineAvailable ? 'IQ-TREE detected. Maximum likelihood and built-in methods are available.' : 'Built-in methods are ready. Install IQ-TREE to enable maximum likelihood.';
  if (inputType !== type) {
    inputText=''; $('file').value=''; $('file-name').textContent='No file selected';
  }
}
$('data-type').addEventListener('change', updateRecipes);
$('demo').addEventListener('click',()=>{
  inputType=$('data-type').value; inputText=examples[inputType];
  $('file').value=''; $('file-name').textContent='Synthetic demonstration';
});
$('file').addEventListener('change',async()=>{
  const file=$('file').files[0];
  if (!file) return;
  const type=$('data-type').value;
  if(file.size>10*1024*1024){$('status').textContent='Use the CLI for files larger than 10 MB.'; return;}
  inputText=await file.text(); inputType=type; $('file-name').textContent=file.name;
});
$('run-form').addEventListener('submit',async event=>{
  event.preventDefault(); $('status').className='';
  if(!inputText || inputType!==$('data-type').value){$('status').textContent='Upload a matching file or load an example first.';return;}
  $('run').disabled=true; $('download').hidden=true; $('tree-wrap').hidden=true;
  $('metrics').hidden=true; $('report-section').hidden=true; $('warnings').replaceChildren();
  $('placeholder').hidden=false; $('status').textContent='Running analysis on this machine…';
  try {
    const response=await fetch('/api/infer',{method:'POST',headers:{'Content-Type':'application/json','X-MethylPhylo-Request':'1'},body:JSON.stringify({text:inputText,data_type:inputType,recipe:$('recipe').value,missing:$('missing').value,sites:$('sites').value,seed:Number($('seed').value),bootstrap:Number($('bootstrap').value)})});
    const data=await response.json(); if(!response.ok) throw new Error(data.error || 'Analysis failed');
    const report=data.report;
    $('placeholder').hidden=true; $('tree-wrap').hidden=false;
    if(svgUrl) URL.revokeObjectURL(svgUrl);
    svgUrl=URL.createObjectURL(new Blob([data.svg],{type:'image/svg+xml'})); $('tree').src=svgUrl;
    $('metrics').replaceChildren(); $('metrics').hidden=false;
    for(const [label,value] of [['Samples',report.data.samples.length],['Retained sites',report.data.sites],['Inference seconds',report.inference_seconds.toFixed(2)]]) {
      const div=document.createElement('div'); div.className='metric';
      const strong=document.createElement('strong'); strong.textContent=value;
      div.append(strong,document.createTextNode(label)); $('metrics').append(div);
    }
    for(const warning of report.warnings){const li=document.createElement('li');li.textContent=warning;$('warnings').append(li);}
    $('report').textContent=JSON.stringify(report,null,2);$('report-section').hidden=false;
    if(zipUrl) URL.revokeObjectURL(zipUrl);
    zipUrl=URL.createObjectURL(new Blob([Uint8Array.from(atob(data.zip),c=>c.charCodeAt(0))],{type:'application/zip'}));
    $('download').href=zipUrl;$('download').hidden=false;
    $('status').textContent=`Completed: ${report.label}. Download this run to retain all inputs and settings.`;
  } catch(error){$('status').className='error';$('status').textContent=error.message;}
  finally{$('run').disabled=false;}
});
fetch('/api/catalog').then(r=>{if(!r.ok) throw new Error();return r.json();}).then(data=>{recipes=data.recipes;engineAvailable=data.iqtree_available;updateRecipes();}).catch(()=>{$('status').textContent='Start the local workbench with methylphylo-research serve.';$('run').disabled=true;});
