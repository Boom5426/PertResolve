const dotCloud=(kind='perturb')=>`<div class="mini-cloud ${kind}" aria-hidden="true"><i></i><i></i><i></i><i></i><i></i><i></i><i></i><i></i><i></i></div>`;

const visuals={
  detection:()=>`
    <div class="diag-vis detection-vis">
      <div class="diag-vis-heading"><span>Same perturbation, two disjoint cell halves</span><span>Compare with reference</span></div>
      <div class="detect-grid">
        <div class="compare-block">
          <div class="compare-stage">
            <div class="cloud-wrap"><span>half A</span>${dotCloud('perturb')}</div>
            <div class="distance-tag"><b>d<sub>self</sub></b><span>split-half noise</span></div>
            <div class="cloud-wrap"><span>half B</span>${dotCloud('perturb alt')}</div>
          </div>
        </div>
        <div class="compare-divider">→</div>
        <div class="compare-block">
          <div class="compare-stage">
            <div class="cloud-wrap"><span>perturbation</span>${dotCloud('perturb')}</div>
            <div class="distance-tag strong"><b>d<sub>null</sub></b><span>effect vs reference</span></div>
            <div class="cloud-wrap"><span>reference</span>${dotCloud('control')}</div>
          </div>
        </div>
      </div>
      <div class="criterion-band">
        <span><b>signal</b> = d<sub>null</sub> − d<sub>self</sub></span>
        <span><b>width</b> = 95% spread of d<sub>self</sub> across splits</span>
        <strong>detectable when signal &gt; width</strong>
      </div>
      <div class="vis-caption">The question is not “is the mean shifted?” but whether the perturbation–reference separation exceeds the perturbation’s own sampling variability.</div>
    </div>`,
  identification:()=>`
    <div class="diag-vis identification-vis">
      <div class="diag-vis-heading"><span>Local competition among biologically relevant candidates</span><span>Detection is not enough</span></div>
      <div class="candidate-scene">
        <div class="candidate-node target"><span>target allele</span>${dotCloud('perturb')}</div>
        <div class="candidate-gap"><b>local separation</b><span>nearest competitor</span></div>
        <div class="candidate-node sibling"><span>sibling allele</span>${dotCloud('sibling')}</div>
        <div class="candidate-node distant"><span>unrelated perturbation</span>${dotCloud('distant')}</div>
      </div>
      <div class="criterion-band">
        <span><b>separation</b> = cross-fitted nearest-competitor distance</span>
        <span><b>spread</b> = its own 95% sampling spread</span>
        <strong>identifiable when separation &gt; spread</strong>
      </div>
      <div class="vis-caption">A perturbation may be far from control yet still overlap its nearest sibling. Fine-grained scoring depends on the closest confusable alternative.</div>
    </div>`,
  reproducibility:()=>`
    <div class="diag-vis reproducibility-vis">
      <div class="diag-vis-heading"><span>Disjoint measurement A</span><span>Disjoint measurement B</span></div>
      <div class="repro-grid">
        <div class="repro-side">
          <div class="profile-chip"><span>A</span><b>variant 1</b></div>
          <div class="profile-chip"><span>B</span><b>variant 2</b></div>
          <div class="profile-chip"><span>C</span><b>variant 3</b></div>
        </div>
        <div class="match-lines"><i></i><i></i><i></i><span>same identity?</span></div>
        <div class="repro-side right">
          <div class="profile-chip"><span>A′</span><b>variant 1</b></div>
          <div class="profile-chip"><span>B′</span><b>variant 2</b></div>
          <div class="profile-chip"><span>C′</span><b>variant 3</b></div>
        </div>
      </div>
      <div class="criterion-band compact">
        <span>query = one split</span><span>reference = the other split</span><strong>measure how well the experiment recovers itself</strong>
      </div>
      <div class="vis-caption">This is an empirical reproducibility reference from the measurement itself. It is informative context for prediction, not a universal model-performance ceiling.</div>
    </div>`,
  ranking:()=>`
    <div class="diag-vis ranking-vis">
      <div class="diag-vis-heading"><span>Observed predictor gap</span><span>Repeated benchmark measurements</span></div>
      <div class="ranking-grid">
        <div class="score-side">
          <div class="score-row"><span>Model A</span><div class="score-track"><i style="--left:58%;--width:25%"></i><b style="--x:72%"></b></div></div>
          <div class="score-row"><span>Model B</span><div class="score-track"><i style="--left:42%;--width:27%"></i><b style="--x:56%"></b></div></div>
          <div class="score-gap"><span>Δ model score</span><strong>is this gap larger than ranking uncertainty?</strong></div>
        </div>
        <div class="rank-outcomes">
          <div><span class="stable-dot"></span><b>A &gt; B</b><small>stable repeats</small></div>
          <div><span class="mixed-dot"></span><b>A ≈ B</b><small>order flips</small></div>
        </div>
      </div>
      <div class="criterion-band compact">
        <span>depends on score gap</span><span>depends on benchmark size</span><strong>ranking resolution is separate from response detectability</strong>
      </div>
      <div class="vis-caption">A benchmark can measure perturbation responses yet still be too noisy or too small to reliably order two predictors whose scores differ only slightly.</div>
    </div>`
};

const panels={
  detection:{
    number:'01 / Detection',
    title:'Can the perturbation be distinguished from reference?',
    copy:'Two disjoint halves of the same perturbation define its split-half variability (d_self), while perturbation versus reference defines d_null. Detection asks whether the effect clears that perturbation’s own measurement noise.',
    note:'Detection answers whether there is a reproducibly measurable response. It does not imply that nearby perturbations can be identified.',
    visual:visuals.detection
  },
  identification:{
    number:'02 / Identification',
    title:'Can the response be distinguished from its relevant competitors?',
    copy:'Identification evaluates local separation from the nearest biologically relevant competitor, such as another allele of the same gene, using disjoint cells from those used for detection.',
    note:'A benchmark can support broad perturbation detection while remaining unresolved for fine-grained allele identification.',
    visual:visuals.identification
  },
  reproducibility:{
    number:'03 / Split-half reference',
    title:'Does the measured response reproduce across disjoint samples?',
    copy:'Independent cell splits provide an empirical reference for how consistently the experiment recovers the same perturbation-level structure.',
    note:'This reference characterizes the measurement. It is not a universal model-performance ceiling.',
    visual:visuals.reproducibility
  },
  ranking:{
    number:'04 / Model-ranking resolution',
    title:'Which predictor differences can the benchmark reliably order?',
    copy:'Controlled predictor differences are compared against benchmark variability to ask which score gaps remain reproducibly rankable.',
    note:'Ranking resolution depends on the predictor gap, measured signal, candidate geometry, sampling depth and benchmark size.',
    visual:visuals.ranking
  }
};

const diagnosticStyle=document.createElement('style');
diagnosticStyle.id='diagnostic-visual-styles';
diagnosticStyle.textContent=`
.evidence-graphic{height:auto;min-height:250px;margin:26px 0 18px;border:0!important;background:none!important}
.evidence-graphic:before,.evidence-graphic:after{display:none!important}
.evidence-graphic>.axis-x,.evidence-graphic>.axis-y{display:none!important}
.diag-vis{border-top:1px solid #1e1e1e;border-bottom:1px solid #d8d4cc;background:#fff;padding:13px 0 0;color:#343434}
.diag-vis-heading{display:flex;justify-content:space-between;gap:20px;padding:0 3px 10px;font-size:8px;font-weight:700;letter-spacing:.12em;text-transform:uppercase;color:#6f6d69}
.detect-grid{display:grid;grid-template-columns:1fr 24px 1fr;align-items:center;border-top:1px solid #e2ded6;padding:24px 10px 20px}
.compare-stage{display:grid;grid-template-columns:1fr 88px 1fr;align-items:center;gap:7px}
.cloud-wrap{text-align:center}.cloud-wrap>span,.candidate-node>span{display:block;margin-bottom:8px;font-size:8px;color:#6f6d69;text-transform:uppercase;letter-spacing:.08em}
.mini-cloud{position:relative;width:82px;height:68px;margin:auto;border:1px solid #dedad2;border-radius:50%;background:#faf9f6}
.mini-cloud i{position:absolute;width:7px;height:7px;border-radius:50%;background:#b43a32;opacity:.78}
.mini-cloud i:nth-child(1){left:22%;top:22%}.mini-cloud i:nth-child(2){left:48%;top:17%}.mini-cloud i:nth-child(3){left:66%;top:31%}.mini-cloud i:nth-child(4){left:30%;top:48%}.mini-cloud i:nth-child(5){left:54%;top:45%}.mini-cloud i:nth-child(6){left:72%;top:57%}.mini-cloud i:nth-child(7){left:18%;top:68%}.mini-cloud i:nth-child(8){left:46%;top:73%}.mini-cloud i:nth-child(9){left:62%;top:67%}
.mini-cloud.alt i{transform:translate(2px,-1px)}.mini-cloud.control i{background:#626a70}.mini-cloud.sibling i{background:#d1934b}.mini-cloud.distant i{background:#7b8790}
.distance-tag{text-align:center}.distance-tag:before{content:'';display:block;height:1px;background:#9c9991;margin-bottom:6px}.distance-tag b{display:block;font:400 17px/1 Georgia,serif;color:#1e1e1e}.distance-tag span{display:block;margin-top:4px;font-size:7px;line-height:1.3;color:#6f6d69}.distance-tag.strong:before{height:2px;background:#b43a32}.distance-tag.strong b{color:#b43a32}
.compare-divider{text-align:center;color:#b43a32;font-size:18px}
.criterion-band{display:grid;grid-template-columns:1fr 1.35fr 1.2fr;gap:0;border-top:1px solid #d8d4cc;background:#f6f4ef}.criterion-band>*{padding:11px 12px;border-right:1px solid #d8d4cc;font-size:8px;line-height:1.45}.criterion-band>*:last-child{border-right:0}.criterion-band strong{color:#8f2e28}.criterion-band.compact{grid-template-columns:1fr 1fr 1.65fr}
.vis-caption{padding:10px 3px 12px;font-size:9px;line-height:1.55;color:#6f6d69}
.candidate-scene{position:relative;display:grid;grid-template-columns:1fr 115px 1fr 1fr;gap:13px;align-items:end;padding:25px 8px 22px;border-top:1px solid #e2ded6}.candidate-node{text-align:center}.candidate-node.distant{opacity:.64;transform:translateX(5px)}.candidate-gap{text-align:center;align-self:center}.candidate-gap:before{content:'';display:block;height:2px;background:#b43a32;margin-bottom:7px}.candidate-gap b{display:block;font-size:9px;color:#8f2e28}.candidate-gap span{display:block;margin-top:3px;font-size:7px;color:#6f6d69}
.repro-grid{display:grid;grid-template-columns:1fr 90px 1fr;gap:18px;align-items:center;padding:20px 16px;border-top:1px solid #e2ded6}.repro-side{display:grid;gap:8px}.profile-chip{display:grid;grid-template-columns:28px 1fr;align-items:center;border:1px solid #d8d4cc;background:#faf9f6}.profile-chip span{display:grid;place-items:center;height:28px;border-right:1px solid #d8d4cc;font:400 15px Georgia,serif;color:#b43a32}.profile-chip b{padding-left:10px;font-size:9px}.repro-side.right .profile-chip span{color:#626a70}.match-lines{display:grid;gap:14px;align-items:center}.match-lines i{display:block;height:1px;background:#9aa0a5}.match-lines span{font-size:7px;text-align:center;color:#6f6d69;text-transform:uppercase;letter-spacing:.08em}
.ranking-grid{display:grid;grid-template-columns:1.45fr .55fr;gap:28px;padding:22px 14px;border-top:1px solid #e2ded6;align-items:center}.score-side{display:grid;gap:13px}.score-row{display:grid;grid-template-columns:58px 1fr;gap:10px;align-items:center;font-size:9px}.score-track{position:relative;height:16px;background:#efede8}.score-track i{position:absolute;left:var(--left);width:var(--width);top:5px;height:6px;background:#c8c4bc}.score-track b{position:absolute;left:var(--x);top:2px;width:2px;height:12px;background:#b43a32}.score-gap{display:flex;justify-content:space-between;gap:12px;border-top:1px solid #d8d4cc;padding-top:9px;font-size:8px;color:#6f6d69}.score-gap strong{color:#343434}.rank-outcomes{display:grid;gap:10px}.rank-outcomes>div{display:grid;grid-template-columns:12px 1fr;column-gap:6px;align-items:center;border-bottom:1px solid #e2ded6;padding-bottom:8px}.rank-outcomes b{font-size:9px}.rank-outcomes small{grid-column:2;font-size:7px;color:#6f6d69}.stable-dot,.mixed-dot{width:8px;height:8px;border-radius:50%;background:#b43a32}.mixed-dot{background:linear-gradient(90deg,#b43a32 50%,#626a70 50%)}
@media(max-width:650px){.detect-grid{grid-template-columns:1fr}.compare-divider{transform:rotate(90deg);margin:4px 0}.compare-stage{grid-template-columns:1fr 70px 1fr}.criterion-band,.criterion-band.compact{grid-template-columns:1fr}.criterion-band>*{border-right:0;border-bottom:1px solid #d8d4cc}.criterion-band>*:last-child{border-bottom:0}.candidate-scene{grid-template-columns:1fr 75px 1fr}.candidate-node.distant{display:none}.repro-grid{grid-template-columns:1fr 45px 1fr;gap:8px;padding-left:4px;padding-right:4px}.ranking-grid{grid-template-columns:1fr}.mini-cloud{width:68px;height:58px}.diag-vis-heading{font-size:7px}}
`;
document.head.appendChild(diagnosticStyle);

function renderPanel(key){
  const panel=panels[key];
  if(!panel)return;
  document.querySelector('#panel-number').textContent=panel.number;
  document.querySelector('#panel-title').textContent=panel.title;
  document.querySelector('#panel-copy').textContent=panel.copy;
  document.querySelector('#panel-note').textContent=panel.note;
  const visual=document.querySelector('.evidence-graphic');
  if(visual)visual.innerHTML=panel.visual();
}

document.querySelectorAll('.diagnostic').forEach(button=>{
  button.addEventListener('click',()=>{
    document.querySelectorAll('.diagnostic').forEach(item=>{
      item.classList.toggle('active',item===button);
      item.setAttribute('aria-selected',String(item===button));
    });
    renderPanel(button.dataset.panel);
  });
});
renderPanel('detection');

const menuButton=document.querySelector('.menu-toggle');
const nav=document.querySelector('#site-nav');
menuButton?.addEventListener('click',()=>{
  const open=nav.classList.toggle('is-open');
  menuButton.setAttribute('aria-expanded',String(open));
});
nav?.querySelectorAll('a').forEach(link=>link.addEventListener('click',()=>{
  nav.classList.remove('is-open');
  menuButton?.setAttribute('aria-expanded','false');
}));

document.querySelectorAll('.copy-button').forEach(button=>button.addEventListener('click',async()=>{
  const target=document.querySelector(`#${button.dataset.copyTarget}`);
  if(!target)return;
  try{
    await navigator.clipboard.writeText(target.innerText);
    button.textContent='Copied';
    setTimeout(()=>button.textContent='Copy',1600);
  }catch{
    button.textContent='Select code';
  }
}));