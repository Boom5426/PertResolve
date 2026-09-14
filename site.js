const panels = {
  detection: {
    number: '01 / detection',
    title: 'Is there a detectable response?',
    copy: 'PertResolve compares each perturbation with a reference population while accounting for split-half sampling noise.',
    note: 'A detectable signal alone does not establish that nearby perturbations, or competing models, can be distinguished.'
  },
  identification: {
    number: '02 / identification',
    title: 'Can nearby responses be separated?',
    copy: 'For each perturbation, the diagnostic asks whether its measured response is distinct from its nearest competitor in the relevant candidate pool.',
    note: 'A broad cross-gene comparison can look easy while a biologically relevant same-gene allele comparison remains unresolved.'
  },
  reproducibility: {
    number: '03 / split-half reference',
    title: 'Does the response reproduce?',
    copy: 'Disjoint subsets of cells provide an empirical reference for the reproducibility of each measured perturbation response.',
    note: 'This is a reference for interpretation, not a hard performance ceiling and not evidence of model-ranking resolution by itself.'
  },
  ranking: {
    number: '04 / model-ranking resolution',
    title: 'Which model gaps can the data order?',
    copy: 'Known predictor differences are tested against the measurement to determine which score gaps can be reliably ranked.',
    note: 'The answer depends on the candidate pool, sampling depth and the actual predictor differences under comparison.'
  }
};

const diagnostics = document.querySelectorAll('.diagnostic');
diagnostics.forEach((item) => {
  item.addEventListener('click', () => {
    const panel = panels[item.dataset.panel];
    diagnostics.forEach((entry) => entry.classList.toggle('active', entry === item));
    document.querySelector('#panel-number').textContent = panel.number;
    document.querySelector('#panel-title').textContent = panel.title;
    document.querySelector('#panel-copy').textContent = panel.copy;
    document.querySelector('#panel-note').textContent = panel.note;
  });
});

const menuButton = document.querySelector('.menu-toggle');
const navigation = document.querySelector('#site-nav');
menuButton.addEventListener('click', () => {
  const isOpen = navigation.classList.toggle('is-open');
  menuButton.setAttribute('aria-expanded', String(isOpen));
});
navigation.querySelectorAll('a').forEach((link) => link.addEventListener('click', () => {
  navigation.classList.remove('is-open');
  menuButton.setAttribute('aria-expanded', 'false');
}));

document.querySelectorAll('.copy-button').forEach((button) => {
  button.addEventListener('click', async () => {
    const code = document.querySelector(`#${button.dataset.copyTarget}`).innerText;
    try {
      await navigator.clipboard.writeText(code);
      button.textContent = 'Copied';
      setTimeout(() => { button.textContent = 'Copy'; }, 1800);
    } catch {
      button.textContent = 'Select code';
    }
  });
});
