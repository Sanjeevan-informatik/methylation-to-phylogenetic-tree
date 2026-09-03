(() => {
  "use strict";

  const GROUP_COLORS = ["#6ee7d8", "#ff8f78", "#7db7ff", "#f6ca68", "#c6a5ff", "#99df7b"];
  const METRIC_HELP = {
    euclidean: "Sensitive to the magnitude of methylation differences.",
    manhattan: "Robust summary of absolute beta-value differences.",
    correlation: "Compares methylation patterns rather than absolute levels.",
  };

  const els = {
    demoButton: document.querySelector("#demo-button"),
    fileInput: document.querySelector("#file-input"),
    fileStatus: document.querySelector("#file-status"),
    resetButton: document.querySelector("#reset-button"),
    metricSelect: document.querySelector("#metric-select"),
    metricHelp: document.querySelector("#metric-help"),
    lociRange: document.querySelector("#loci-range"),
    lociOutput: document.querySelector("#loci-output"),
    rangeMax: document.querySelector("#range-max"),
    runButton: document.querySelector("#run-button"),
    errorMessage: document.querySelector("#error-message"),
    newickButton: document.querySelector("#newick-button"),
    reportButton: document.querySelector("#report-button"),
    metricSamples: document.querySelector("#metric-samples"),
    metricGroups: document.querySelector("#metric-groups"),
    metricLoci: document.querySelector("#metric-loci"),
    metricSelected: document.querySelector("#metric-selected"),
    metricMissing: document.querySelector("#metric-missing"),
    metricHeight: document.querySelector("#metric-height"),
    treeSubtitle: document.querySelector("#tree-subtitle"),
    pcaSubtitle: document.querySelector("#pca-subtitle"),
    legend: document.querySelector("#legend"),
    treeChart: document.querySelector("#tree-chart"),
    pcaChart: document.querySelector("#pca-chart"),
    heatmapChart: document.querySelector("#heatmap-chart"),
    qcChecks: document.querySelector("#qc-checks"),
    distanceTable: document.querySelector("#distance-table"),
    modelSelect: document.querySelector("#model-select"),
    modelDistance: document.querySelector("#model-distance"),
    modelDistanceLabel: document.querySelector("#model-distance-label"),
    modelContext: document.querySelector("#model-context"),
    modelResultTitle: document.querySelector("#model-result-title"),
    modelResultCaption: document.querySelector("#model-result-caption"),
    probSame: document.querySelector("#prob-same"),
    probChange: document.querySelector("#prob-change"),
    probSpecific: document.querySelector("#prob-specific"),
    probChangeLabel: document.querySelector("#prob-change-label"),
    probSameBar: document.querySelector("#prob-same-bar"),
    probChangeBar: document.querySelector("#prob-change-bar"),
    specificProbability: document.querySelector("#specific-probability"),
    modelFormula: document.querySelector("#model-formula"),
    probabilityCheck: document.querySelector("#probability-check"),
  };

  let source = makeDemoData();
  let sourceLabel = "Synthetic AML cohort";
  let latest = null;

  function makeDemoData() {
    const samples = [
      "Normal_01", "Normal_02", "Normal_03", "Primary_01",
      "Primary_02", "Relapse_01", "Relapse_02", "Relapse_03",
    ];
    const rows = [];
    for (let locus = 0; locus < 72; locus += 1) {
      const baseline = 0.18 + 0.55 * ((Math.sin(locus * 1.71) + 1) / 2);
      const values = samples.map((sample, sampleIndex) => {
        const group = sample.split("_")[0];
        let shift = 0;
        if (group === "Primary") {
          shift = locus % 3 === 0 ? 0.23 : locus % 4 === 0 ? -0.17 : 0.045;
        } else if (group === "Relapse") {
          shift = locus % 3 === 0 ? 0.31 : locus % 4 === 0 ? -0.24 : locus % 5 === 0 ? 0.17 : 0.075;
        }
        const noise = 0.018 * Math.sin((locus + 2) * (sampleIndex + 1) * 0.83);
        return Math.max(0.015, Math.min(0.985, baseline + shift + noise));
      });
      rows.push({ id: `cg${String(10200000 + locus * 137).padStart(8, "0")}`, values });
    }
    [[5, 1], [12, 5], [18, 3], [28, 0], [37, 6], [43, 2], [58, 7], [66, 4]].forEach(([row, col]) => {
      rows[row].values[col] = null;
    });
    return { samples, rows };
  }

  function parseMatrix(text) {
    const clean = text.replace(/^\uFEFF/, "").trim();
    if (!clean) throw new Error("The selected file is empty.");
    const lines = clean.split(/\r?\n/).filter((line) => line.trim());
    if (lines.length < 6) throw new Error("The matrix needs a header and at least five CpG rows.");
    const first = lines[0];
    const delimiter = (first.match(/\t/g) || []).length >= (first.match(/,/g) || []).length ? "\t" : ",";
    const header = parseDelimitedLine(first, delimiter).map((item) => item.trim());
    if (header.length < 4) throw new Error("Provide at least three sample columns after the CpG ID column.");
    const samples = header.slice(1);
    if (new Set(samples).size !== samples.length || samples.some((name) => !name)) {
      throw new Error("Every sample column needs a unique, non-empty name.");
    }
    if (samples.length > 120) throw new Error("For an interactive analysis, use no more than 120 samples.");

    const rows = lines.slice(1).map((line, lineIndex) => {
      const cells = parseDelimitedLine(line, delimiter);
      if (cells.length !== header.length) {
        throw new Error(`Row ${lineIndex + 2} has ${cells.length} columns; expected ${header.length}.`);
      }
      const id = cells[0].trim();
      if (!id) throw new Error(`Row ${lineIndex + 2} has no CpG identifier.`);
      const values = cells.slice(1).map((cell) => {
        const value = cell.trim();
        if (!value || /^(na|nan|null|\.)$/i.test(value)) return null;
        const numeric = Number(value);
        if (!Number.isFinite(numeric)) throw new Error(`Non-numeric beta value “${value}” at ${id}.`);
        if (numeric < 0 || numeric > 1) throw new Error(`Beta value ${numeric} at ${id} is outside 0–1.`);
        return numeric;
      });
      return { id, values };
    });
    return { samples, rows };
  }

  function parseDelimitedLine(line, delimiter) {
    const cells = [];
    let current = "";
    let quoted = false;
    for (let index = 0; index < line.length; index += 1) {
      const character = line[index];
      if (character === '"') {
        if (quoted && line[index + 1] === '"') {
          current += '"';
          index += 1;
        } else {
          quoted = !quoted;
        }
      } else if (character === delimiter && !quoted) {
        cells.push(current);
        current = "";
      } else {
        current += character;
      }
    }
    cells.push(current);
    return cells;
  }

  function median(values) {
    const sorted = [...values].sort((a, b) => a - b);
    const middle = Math.floor(sorted.length / 2);
    return sorted.length % 2 ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2;
  }

  function variance(values) {
    const mean = values.reduce((sum, value) => sum + value, 0) / values.length;
    return values.reduce((sum, value) => sum + (value - mean) ** 2, 0) / values.length;
  }

  function preprocess(matrix, selectedCount) {
    let missing = 0;
    let dropped = 0;
    const cleaned = [];
    matrix.rows.forEach((row) => {
      const observed = row.values.filter((value) => value !== null);
      missing += row.values.length - observed.length;
      if (!observed.length) {
        dropped += 1;
        return;
      }
      const fill = median(observed);
      const values = row.values.map((value) => (value === null ? fill : value));
      cleaned.push({ id: row.id, values, variance: variance(values) });
    });
    if (cleaned.length < 5) throw new Error("At least five CpG loci must contain observed values.");
    cleaned.sort((a, b) => b.variance - a.variance);
    const count = Math.max(5, Math.min(selectedCount, cleaned.length));
    const selected = cleaned.slice(0, count);
    return {
      selected,
      all: cleaned,
      missing,
      dropped,
      missingRate: missing / (matrix.rows.length * matrix.samples.length),
      lowVariance: cleaned.filter((row) => row.variance < 0.0001).length,
    };
  }

  function pairwiseDistances(selected, sampleCount, metric) {
    const distances = Array.from({ length: sampleCount }, () => Array(sampleCount).fill(0));
    for (let left = 0; left < sampleCount; left += 1) {
      for (let right = left + 1; right < sampleCount; right += 1) {
        const a = selected.map((row) => row.values[left]);
        const b = selected.map((row) => row.values[right]);
        let distance;
        if (metric === "manhattan") {
          distance = a.reduce((sum, value, index) => sum + Math.abs(value - b[index]), 0) / a.length;
        } else if (metric === "correlation") {
          const meanA = a.reduce((sum, value) => sum + value, 0) / a.length;
          const meanB = b.reduce((sum, value) => sum + value, 0) / b.length;
          let numerator = 0;
          let sumA = 0;
          let sumB = 0;
          a.forEach((value, index) => {
            const da = value - meanA;
            const db = b[index] - meanB;
            numerator += da * db;
            sumA += da ** 2;
            sumB += db ** 2;
          });
          distance = sumA && sumB ? 1 - numerator / Math.sqrt(sumA * sumB) : 1;
        } else {
          distance = Math.sqrt(a.reduce((sum, value, index) => sum + (value - b[index]) ** 2, 0) / a.length);
        }
        distances[left][right] = distance;
        distances[right][left] = distance;
      }
    }
    return distances;
  }

  function upgma(samples, distances) {
    let nextId = samples.length;
    let clusters = samples.map((name, index) => ({
      id: index,
      members: [index],
      node: { name, height: 0, sampleIndex: index, children: null },
    }));

    const clusterDistance = (a, b) => {
      let total = 0;
      let count = 0;
      a.members.forEach((left) => b.members.forEach((right) => {
        total += distances[left][right];
        count += 1;
      }));
      return total / count;
    };

    while (clusters.length > 1) {
      let best = { left: 0, right: 1, distance: clusterDistance(clusters[0], clusters[1]) };
      for (let left = 0; left < clusters.length; left += 1) {
        for (let right = left + 1; right < clusters.length; right += 1) {
          const distance = clusterDistance(clusters[left], clusters[right]);
          if (distance < best.distance) best = { left, right, distance };
        }
      }
      const a = clusters[best.left];
      const b = clusters[best.right];
      const merged = {
        id: nextId,
        members: [...a.members, ...b.members],
        node: {
          name: null,
          height: best.distance / 2,
          children: [a.node, b.node],
        },
      };
      nextId += 1;
      clusters = clusters.filter((_, index) => index !== best.left && index !== best.right);
      clusters.push(merged);
    }
    return clusters[0].node;
  }

  function newick(node, parentHeight = null) {
    const branch = parentHeight === null ? "" : `:${Math.max(0, parentHeight - node.height).toFixed(6)}`;
    if (!node.children) return `${safeNewickLabel(node.name)}${branch}`;
    return `(${node.children.map((child) => newick(child, node.height)).join(",")})${branch}`;
  }

  function safeNewickLabel(label) {
    return /^[A-Za-z0-9_.-]+$/.test(label) ? label : `'${label.replaceAll("'", "''")}'`;
  }

  function inferGroups(samples) {
    const groups = [];
    const map = new Map();
    samples.forEach((sample) => {
      const match = sample.match(/^(.+?)[_\-. ]?\d+$/);
      const group = match && match[1].length > 1 ? match[1] : "Samples";
      if (!map.has(group)) {
        map.set(group, { name: group, color: GROUP_COLORS[groups.length % GROUP_COLORS.length], samples: [] });
        groups.push(map.get(group));
      }
      map.get(group).samples.push(sample);
    });
    if (groups.length === samples.length || groups.every((group) => group.samples.length === 1)) {
      return [{ name: "Samples", color: GROUP_COLORS[0], samples: [...samples] }];
    }
    return groups;
  }

  function colorForSample(sample, groups) {
    return groups.find((group) => group.samples.includes(sample))?.color || GROUP_COLORS[0];
  }

  function runAnalysis() {
    clearError();
    try {
      const metric = els.metricSelect.value;
      const processed = preprocess(source, Number(els.lociRange.value));
      const distances = pairwiseDistances(processed.selected, source.samples.length, metric);
      const root = upgma(source.samples, distances);
      const groups = inferGroups(source.samples);
      const pca = calculatePca(processed.selected, source.samples.length);
      latest = { metric, processed, distances, root, groups, pca, newick: `${newick(root)};` };
      renderResults();
    } catch (error) {
      showError(error instanceof Error ? error.message : "The analysis could not be completed.");
    }
  }

  function renderResults() {
    const { metric, processed, distances, root, groups, pca } = latest;
    els.metricSamples.textContent = String(source.samples.length);
    els.metricGroups.textContent = groups.length === 1 ? "unlabeled samples" : `${groups.length} inferred groups`;
    els.metricLoci.textContent = String(source.rows.length);
    els.metricSelected.textContent = `${processed.selected.length} selected`;
    els.metricMissing.textContent = `${(processed.missingRate * 100).toFixed(1)}%`;
    els.metricHeight.textContent = root.height.toFixed(3);
    els.treeSubtitle.textContent = `${metricLabel(metric)} · ${processed.selected.length} high-variance loci`;
    els.pcaSubtitle.textContent = `PC1 ${pca.explained[0].toFixed(1)}% · PC2 ${pca.explained[1].toFixed(1)}% of sample-space variance`;
    els.legend.innerHTML = groups.map((group) => `<span><i style="background:${group.color}"></i>${escapeHtml(group.name)}</span>`).join("");
    renderTree(root, groups);
    renderPca(pca, groups);
    renderHeatmap(processed.selected, groups);
    renderQc(processed);
    renderDistanceTable(distances);
  }

  function metricLabel(metric) {
    return { euclidean: "Euclidean distance", manhattan: "Manhattan distance", correlation: "Correlation distance" }[metric];
  }

  function transitionProbabilities(model, distance) {
    if (!Number.isFinite(distance) || distance < 0) {
      throw new Error("Model distance must be a non-negative number.");
    }
    if (model === "binary") {
      const decay = Math.exp(-2 * distance);
      const same = 0.5 + 0.5 * decay;
      return { model, distance, same, change: 0.5 - 0.5 * decay, specific: null };
    }
    const decay = Math.exp((-4 * distance) / 3);
    const same = 0.25 + 0.75 * decay;
    const specific = 0.25 - 0.25 * decay;
    return { model: "jc69", distance, same, change: 3 * specific, specific };
  }

  function renderModelProbability() {
    const model = els.modelSelect.value;
    const distance = Number(els.modelDistance.value);
    try {
      const probability = transitionProbabilities(model, distance);
      els.modelDistance.removeAttribute("aria-invalid");
      els.modelDistance.title = "";
      const samePercent = 100 * probability.same;
      const changePercent = 100 * probability.change;
      els.probSame.textContent = `${samePercent.toFixed(2)}%`;
      els.probChange.textContent = `${changePercent.toFixed(2)}%`;
      els.probSameBar.style.width = `${Math.min(100, samePercent)}%`;
      els.probChangeBar.style.width = `${Math.min(100, changePercent)}%`;
      els.modelResultCaption.textContent = `At model distance d = ${distance.toFixed(3)}`;
      els.probabilityCheck.textContent = `Probability check: ${samePercent.toFixed(2)}% + ${changePercent.toFixed(2)}% = ${(samePercent + changePercent).toFixed(2)}%`;

      if (model === "binary") {
        els.modelDistanceLabel.textContent = "Distance d · transition rate × time";
        els.modelContext.textContent = "The symmetric two-state model represents methylated ↔ unmethylated switching with equal forward and reverse rates.";
        els.modelResultTitle.textContent = "Two-state transition probabilities";
        els.probChangeLabel.textContent = "Switched state";
        els.specificProbability.hidden = true;
        els.modelFormula.innerHTML = "P(same) = ½ + ½e<sup>−2d</sup>";
      } else {
        els.modelDistanceLabel.textContent = "Distance d · substitutions per site";
        els.modelContext.textContent = "JC69 models A, C, G and T with equal frequencies and equal substitution rates.";
        els.modelResultTitle.textContent = "JC69 transition probabilities";
        els.probChangeLabel.textContent = "Any different nucleotide";
        els.specificProbability.hidden = false;
        els.probSpecific.textContent = `${(100 * probability.specific).toFixed(2)}%`;
        els.modelFormula.innerHTML = "P(same) = ¼ + ¾e<sup>−4d/3</sup>";
      }
    } catch (error) {
      els.modelDistance.setAttribute("aria-invalid", "true");
      els.modelDistance.title = error instanceof Error ? error.message : "Invalid model distance";
      els.probSame.textContent = "—";
      els.probChange.textContent = "—";
      els.probSpecific.textContent = "—";
      els.probSameBar.style.width = "0";
      els.probChangeBar.style.width = "0";
      els.probabilityCheck.textContent = "Enter a non-negative model distance.";
    }
  }

  function renderTree(root, groups) {
    const width = 980;
    const leaves = [];
    const collect = (node) => node.children ? node.children.forEach(collect) : leaves.push(node);
    collect(root);
    const height = Math.max(440, leaves.length * 50 + 70);
    const left = 70;
    const right = 170;
    const top = 44;
    const plotWidth = width - left - right;
    const rootHeight = Math.max(root.height, 0.000001);
    const yStep = (height - top * 2) / Math.max(1, leaves.length - 1);
    leaves.forEach((leaf, index) => { leaf._y = top + index * yStep; });
    const position = (node) => {
      if (node.children) {
        node.children.forEach(position);
        node._y = (node.children[0]._y + node.children[1]._y) / 2;
      }
      node._x = left + ((rootHeight - node.height) / rootHeight) * plotWidth;
    };
    position(root);
    const links = [];
    const nodes = [];
    const walk = (node) => {
      if (node.children) {
        const [a, b] = node.children;
        links.push(`<path class="tree-link" d="M${node._x},${a._y}V${b._y}M${node._x},${a._y}H${a._x}M${node._x},${b._y}H${b._x}"/>`);
        node.children.forEach(walk);
      } else {
        const color = colorForSample(node.name, groups);
        nodes.push(`<circle cx="${node._x}" cy="${node._y}" r="4" fill="${color}"/><text class="tree-label" x="${node._x + 12}" y="${node._y + 4}">${escapeHtml(node.name)}</text>`);
      }
    };
    walk(root);
    const guides = [0, 0.25, 0.5, 0.75, 1].map((fraction) => {
      const x = left + fraction * plotWidth;
      const distance = rootHeight * (1 - fraction);
      return `<line class="tree-guide" x1="${x}" x2="${x}" y1="24" y2="${height - 24}"/><text class="tree-distance" x="${x}" y="18" text-anchor="middle">${distance.toFixed(3)}</text>`;
    }).join("");
    els.treeChart.innerHTML = `<svg viewBox="0 0 ${width} ${height}" role="img" aria-labelledby="tree-title tree-desc"><title id="tree-title">UPGMA sample tree</title><desc id="tree-desc">A dendrogram of ${leaves.length} samples based on ${escapeHtml(metricLabel(latest.metric))}.</desc>${guides}${links.join("")}${nodes.join("")}</svg>`;
  }

  function calculatePca(selected, sampleCount) {
    const samples = Array.from({ length: sampleCount }, (_, sampleIndex) => selected.map((row) => row.values[sampleIndex]));
    for (let feature = 0; feature < selected.length; feature += 1) {
      const mean = samples.reduce((sum, sample) => sum + sample[feature], 0) / sampleCount;
      samples.forEach((sample) => { sample[feature] -= mean; });
    }
    const gram = Array.from({ length: sampleCount }, () => Array(sampleCount).fill(0));
    for (let left = 0; left < sampleCount; left += 1) {
      for (let right = left; right < sampleCount; right += 1) {
        const value = samples[left].reduce((sum, feature, index) => sum + feature * samples[right][index], 0) / Math.max(1, selected.length - 1);
        gram[left][right] = value;
        gram[right][left] = value;
      }
    }
    const first = principalEigen(gram, 0);
    const deflated = gram.map((row, i) => row.map((value, j) => value - first.value * first.vector[i] * first.vector[j]));
    const second = principalEigen(deflated, 1);
    const total = gram.reduce((sum, row, index) => sum + row[index], 0) || 1;
    return {
      points: source.samples.map((name, index) => ({
        name,
        x: first.vector[index] * Math.sqrt(Math.max(0, first.value)),
        y: second.vector[index] * Math.sqrt(Math.max(0, second.value)),
      })),
      explained: [100 * first.value / total, 100 * second.value / total],
    };
  }

  function principalEigen(matrix, seedOffset) {
    const size = matrix.length;
    let vector = Array.from({ length: size }, (_, index) => Math.sin((index + 1) * (1.37 + seedOffset * 0.41)) + 0.3);
    vector = normalize(vector);
    for (let iteration = 0; iteration < 120; iteration += 1) {
      const next = matrix.map((row) => row.reduce((sum, value, index) => sum + value * vector[index], 0));
      const norm = Math.sqrt(next.reduce((sum, value) => sum + value ** 2, 0));
      if (norm < 1e-12) break;
      vector = next.map((value) => value / norm);
    }
    const product = matrix.map((row) => row.reduce((sum, value, index) => sum + value * vector[index], 0));
    const value = vector.reduce((sum, component, index) => sum + component * product[index], 0);
    return { vector, value: Math.max(0, value) };
  }

  function normalize(vector) {
    const norm = Math.sqrt(vector.reduce((sum, value) => sum + value ** 2, 0)) || 1;
    return vector.map((value) => value / norm);
  }

  function renderPca(pca, groups) {
    const width = 980;
    const height = 460;
    const padding = { left: 70, right: 100, top: 35, bottom: 58 };
    const xValues = pca.points.map((point) => point.x);
    const yValues = pca.points.map((point) => point.y);
    const xPad = (Math.max(...xValues) - Math.min(...xValues) || 1) * 0.15;
    const yPad = (Math.max(...yValues) - Math.min(...yValues) || 1) * 0.18;
    const xMin = Math.min(...xValues) - xPad;
    const xMax = Math.max(...xValues) + xPad;
    const yMin = Math.min(...yValues) - yPad;
    const yMax = Math.max(...yValues) + yPad;
    const scaleX = (value) => padding.left + ((value - xMin) / (xMax - xMin)) * (width - padding.left - padding.right);
    const scaleY = (value) => height - padding.bottom - ((value - yMin) / (yMax - yMin)) * (height - padding.top - padding.bottom);
    const grid = [0, 0.25, 0.5, 0.75, 1].map((fraction) => {
      const x = padding.left + fraction * (width - padding.left - padding.right);
      const y = padding.top + fraction * (height - padding.top - padding.bottom);
      return `<line class="pca-grid" x1="${x}" x2="${x}" y1="${padding.top}" y2="${height - padding.bottom}"/><line class="pca-grid" x1="${padding.left}" x2="${width - padding.right}" y1="${y}" y2="${y}"/>`;
    }).join("");
    const points = pca.points.map((point) => {
      const x = scaleX(point.x);
      const y = scaleY(point.y);
      const color = colorForSample(point.name, groups);
      return `<circle cx="${x}" cy="${y}" r="7" fill="${color}" fill-opacity="0.9" stroke="#07131f" stroke-width="2"><title>${escapeHtml(point.name)}</title></circle><text class="pca-sample" x="${x + 11}" y="${y + 4}">${escapeHtml(point.name)}</text>`;
    }).join("");
    els.pcaChart.innerHTML = `<svg viewBox="0 0 ${width} ${height}" role="img" aria-labelledby="pca-title"><title id="pca-title">PCA projection of methylation samples</title>${grid}<line class="pca-axis" x1="${padding.left}" x2="${width - padding.right}" y1="${height - padding.bottom}" y2="${height - padding.bottom}"/><line class="pca-axis" x1="${padding.left}" x2="${padding.left}" y1="${padding.top}" y2="${height - padding.bottom}"/>${points}<text class="pca-label" x="${(padding.left + width - padding.right) / 2}" y="${height - 18}" text-anchor="middle">PC1 · ${pca.explained[0].toFixed(1)}%</text><text class="pca-label" x="18" y="${height / 2}" text-anchor="middle" transform="rotate(-90 18 ${height / 2})">PC2 · ${pca.explained[1].toFixed(1)}%</text></svg>`;
  }

  function renderHeatmap(selected) {
    const rows = selected.slice(0, 18);
    const columns = source.samples.length + 1;
    const header = `<span></span>${source.samples.map((sample) => `<span class="heat-label sample">${escapeHtml(sample)}</span>`).join("")}`;
    const cells = rows.map((row) => `<span class="heat-label" title="${escapeHtml(row.id)}">${escapeHtml(row.id)}</span>${row.values.map((value) => `<span class="heat-cell" title="β = ${value.toFixed(3)}" style="background:${betaColor(value)}"></span>`).join("")}`).join("");
    els.heatmapChart.innerHTML = `<div class="heatmap-grid" style="grid-template-columns:110px repeat(${columns - 1}, minmax(46px, 1fr))">${header}${cells}</div>`;
  }

  function betaColor(value) {
    if (value < 0.5) {
      const t = value * 2;
      return `rgb(${Math.round(17 + 23 * t)},${Math.round(42 + 156 * t)},${Math.round(56 + 124 * t)})`;
    }
    const t = (value - 0.5) * 2;
    return `rgb(${Math.round(40 + 204 * t)},${Math.round(198 + 6 * t)},${Math.round(180 - 68 * t)})`;
  }

  function renderQc(processed) {
    const checks = [
      ["Matrix shape", `${source.rows.length.toLocaleString()} loci × ${source.samples.length} samples parsed successfully.`],
      ["Beta-value range", "All observed values fall within the expected 0–1 interval."],
      ["Missingness", `${processed.missing.toLocaleString()} values (${(processed.missingRate * 100).toFixed(2)}%) replaced with locus medians.`],
      ["Feature selection", `${processed.selected.length} of ${processed.all.length} usable loci retained by variance${processed.dropped ? `; ${processed.dropped} empty rows dropped` : ""}.`],
    ];
    els.qcChecks.innerHTML = checks.map(([title, detail]) => `<div class="qc-row"><span class="qc-icon">✓</span><span><strong>${title}</strong><small>${detail}</small></span></div>`).join("");
  }

  function renderDistanceTable(distances) {
    const header = `<tr><th>Sample</th>${source.samples.map((sample) => `<th title="${escapeHtml(sample)}">${escapeHtml(shortLabel(sample))}</th>`).join("")}</tr>`;
    const rows = distances.map((row, rowIndex) => `<tr><td title="${escapeHtml(source.samples[rowIndex])}">${escapeHtml(shortLabel(source.samples[rowIndex]))}</td>${row.map((value) => `<td>${value.toFixed(3)}</td>`).join("")}</tr>`).join("");
    els.distanceTable.innerHTML = `<table class="distance-table"><thead>${header}</thead><tbody>${rows}</tbody></table>`;
  }

  function shortLabel(label) {
    return label.length > 12 ? `${label.slice(0, 10)}…` : label;
  }

  function updateInputStatus(label, detail, uploaded) {
    sourceLabel = label;
    els.fileStatus.innerHTML = `<span class="status-dot" aria-hidden="true"></span><span><strong>${escapeHtml(label)}</strong><small>${escapeHtml(detail)}</small></span>`;
    els.demoButton.classList.toggle("active", !uploaded);
    const fileTrigger = document.querySelector(".file-trigger");
    fileTrigger.classList.toggle("active", uploaded);
  }

  function updateRange() {
    const max = Math.max(10, Math.min(500, source.rows.length));
    els.lociRange.max = String(max);
    els.rangeMax.textContent = String(max);
    els.lociRange.value = String(Math.min(40, max));
    els.lociOutput.textContent = els.lociRange.value;
  }

  function showError(message) {
    els.errorMessage.textContent = message;
    els.errorMessage.hidden = false;
  }

  function clearError() {
    els.errorMessage.hidden = true;
    els.errorMessage.textContent = "";
  }

  function escapeHtml(value) {
    return String(value).replace(/[&<>"']/g, (character) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#039;" }[character]));
  }

  function download(filename, content, type) {
    const blob = new Blob([content], { type });
    const url = URL.createObjectURL(blob);
    const anchor = document.createElement("a");
    anchor.href = url;
    anchor.download = filename;
    document.body.appendChild(anchor);
    anchor.click();
    anchor.remove();
    URL.revokeObjectURL(url);
  }

  function exportReport() {
    if (!latest) return;
    let modelProbability;
    try {
      modelProbability = transitionProbabilities(els.modelSelect.value, Number(els.modelDistance.value));
    } catch (error) {
      showError(error instanceof Error ? error.message : "The model distance is invalid.");
      return;
    }
    const report = {
      software: "MethylPhylo Studio",
      generated_at: new Date().toISOString(),
      source: sourceLabel,
      parameters: { distance_metric: latest.metric, selected_loci: latest.processed.selected.length, clustering: "UPGMA" },
      quality: {
        samples: source.samples.length,
        input_loci: source.rows.length,
        usable_loci: latest.processed.all.length,
        missing_values: latest.processed.missing,
        missing_rate: Number(latest.processed.missingRate.toFixed(6)),
        empty_loci_dropped: latest.processed.dropped,
      },
      samples: source.samples,
      newick: latest.newick,
      distance_matrix: latest.distances.map((row) => row.map((value) => Number(value.toFixed(6)))),
      transition_probability: modelProbability,
    };
    download("methylphylo-report.json", JSON.stringify(report, null, 2), "application/json");
  }

  document.querySelectorAll("[role='tab']").forEach((tab) => {
    tab.addEventListener("click", () => {
      document.querySelectorAll("[role='tab']").forEach((candidate) => {
        const selected = candidate === tab;
        candidate.classList.toggle("active", selected);
        candidate.setAttribute("aria-selected", String(selected));
        const panel = document.querySelector(`#${candidate.getAttribute("aria-controls")}`);
        panel.hidden = !selected;
        panel.classList.toggle("active", selected);
      });
    });
  });

  els.lociRange.addEventListener("input", () => { els.lociOutput.textContent = els.lociRange.value; });
  els.metricSelect.addEventListener("change", () => { els.metricHelp.textContent = METRIC_HELP[els.metricSelect.value]; });
  els.modelSelect.addEventListener("change", renderModelProbability);
  els.modelDistance.addEventListener("input", renderModelProbability);
  els.runButton.addEventListener("click", runAnalysis);
  els.demoButton.addEventListener("click", () => {
    source = makeDemoData();
    els.fileInput.value = "";
    updateInputStatus("Synthetic AML cohort", "8 samples · 72 CpG loci", false);
    updateRange();
    runAnalysis();
  });
  els.resetButton.addEventListener("click", () => {
    source = makeDemoData();
    els.fileInput.value = "";
    els.metricSelect.value = "euclidean";
    els.metricHelp.textContent = METRIC_HELP.euclidean;
    updateInputStatus("Synthetic AML cohort", "8 samples · 72 CpG loci", false);
    updateRange();
    runAnalysis();
  });
  els.fileInput.addEventListener("change", async () => {
    clearError();
    const file = els.fileInput.files?.[0];
    if (!file) return;
    if (file.size > 8 * 1024 * 1024) {
      showError("Use a CSV or TSV file smaller than 8 MB for the browser demo.");
      return;
    }
    try {
      source = parseMatrix(await file.text());
      updateInputStatus(file.name, `${source.samples.length} samples · ${source.rows.length.toLocaleString()} CpG loci`, true);
      updateRange();
      runAnalysis();
    } catch (error) {
      showError(error instanceof Error ? error.message : "The file could not be parsed.");
    }
  });
  els.newickButton.addEventListener("click", () => latest && download("methylphylo-tree.nwk", `${latest.newick}\n`, "text/plain"));
  els.reportButton.addEventListener("click", exportReport);

  updateRange();
  runAnalysis();
  renderModelProbability();
})();
