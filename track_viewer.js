document.addEventListener('DOMContentLoaded', () => {
    const controls = {
        minExons: document.getElementById('min-exons'),
        maxExons: document.getElementById('max-exons'),
        minTSS: document.getElementById('min-tss'),
        maxTSS: document.getElementById('max-tss'),
        coding: document.getElementById('coding-status'),
        gene: document.getElementById('gene-search')
    };

    const resultsContainer = document.getElementById('results-container');
    const resultsCount = document.getElementById('results-count');
    const tableWrapper = document.getElementById('table-wrapper');
    const applyLink = document.getElementById('apply-filters-ucsc');
    const resetLink = document.getElementById('reset-filters-ucsc');

    const urlParams = new URLSearchParams(window.location.search);
    const hgsid = urlParams.get('hgsid');
    const ucscHost = 'https://genome.ucsc.edu';
    let tableRendered = false;

    // These variables are expected to be defined in the HTML before this script is loaded
    // const transcripts = [...];
    // const trackName = '...';

    function renderTable(data) {
        if (!data || data.length === 0) {
            resultsContainer.innerHTML = '<p>No transcripts match the current filters.</p>';
            return;
        }
        let tableHTML = `
            <table>
                <thead>
                    <tr>
                        <th>Isoform</th><th>Associated Gene</th><th>Exons</th><th>Length</th>
                        <th>Diff to TSS</th><th>Diff to TTS</th><th>Coding</th><th>Isoform Exp</th>
                    </tr>
                </thead>
                <tbody>`;
        data.forEach(t => {
            tableHTML += `
                <tr>
                    <td>${t.isoform || 'N/A'}</td><td>${t.associated_gene || 'N/A'}</td>
                    <td>${t.exons !== null ? t.exons : 'N/A'}</td><td>${t.length !== null ? t.length : 'N/A'}</td>
                    <td>${t.diff_to_TSS !== null ? t.diff_to_TSS : 'N/A'}</td><td>${t.diff_to_TTS !== null ? t.diff_to_TTS : 'N/A'}</td>
                    <td>${t.coding || 'N/A'}</td><td>${t.iso_exp !== null ? t.iso_exp.toFixed(2) : 'N/A'}</td>
                </tr>`;
        });
        tableHTML += '</tbody></table>';
        resultsContainer.innerHTML = tableHTML;
        tableRendered = true;
    }

    function getFilteredIsoforms() {
        const getVal = (el) => el.value;
        const getInt = (el) => el.value ? parseInt(el.value, 10) : null;
        const filters = {
            minExons: getInt(controls.minExons), maxExons: getInt(controls.maxExons),
            minTSS: getInt(controls.minTSS), maxTSS: getInt(controls.maxTSS),
            coding: getVal(controls.coding), gene: getVal(controls.gene).toLowerCase().trim()
        };
        return transcripts.filter(t => {
            if (filters.minExons !== null && (t.exons === null || t.exons < filters.minExons)) return false;
            if (filters.maxExons !== null && (t.exons === null || t.exons > filters.maxExons)) return false;
            if (filters.minTSS !== null && (t.diff_to_TSS === null || t.diff_to_TSS < filters.minTSS)) return false;
            if (filters.maxTSS !== null && (t.diff_to_TSS === null || t.diff_to_TSS > filters.maxTSS)) return false;
            if (filters.coding !== 'all' && t.coding !== filters.coding) return false;
            if (filters.gene && (!t.associated_gene || !t.associated_gene.toLowerCase().includes(filters.gene))) return false;
            return true;
        });
    }
    
    function updateUCSCLink() {
        if (!hgsid) return;
        const filteredData = getFilteredIsoforms();
        const isoform_ids = filteredData.map(t => t.isoform);
        const ucscUrl = new URL(`${ucscHost}/cgi-bin/hgTracks`);
        ucscUrl.searchParams.set('hgsid', hgsid);

        if (isoform_ids.length === 0) {
            ucscUrl.searchParams.set(trackName, 'hide');
        } else {
            const filterText = isoform_ids.join(' ');
            if (filterText.length > 4000) { // Increased limit for modern browsers
                applyLink.classList.add('disabled');
                applyLink.href = '#';
                resultsCount.innerHTML = `Preview: Showing ${filteredData.length} of ${transcripts.length} transcripts. <strong style="color: red;">(Filter too broad to apply to UCSC view)</strong>`;
                return;
            }
            ucscUrl.searchParams.set(`hg_t_${trackName}_filter`, 'on');
            ucscUrl.searchParams.set(`hg_t_${trackName}_filter_text`, filterText);
        }
        applyLink.classList.remove('disabled');
        applyLink.href = ucscUrl.toString();
    }

    function previewFilters() {
        tableWrapper.style.display = 'block';
        const filteredData = getFilteredIsoforms();
        resultsCount.textContent = `Preview: Showing ${filteredData.length} of ${transcripts.length} transcripts.`;
        renderTable(filteredData);
        updateUCSCLink();
    }
    
    function toggleTable() {
        const isHidden = tableWrapper.style.display === 'none';
        tableWrapper.style.display = isHidden ? 'block' : 'none';
        if (isHidden && !tableRendered) {
            renderTable(transcripts);
        }
    }

    // --- Attach Event Listeners ---
    document.getElementById('preview-filters').addEventListener('click', previewFilters);
    document.getElementById('toggle-table').addEventListener('click', toggleTable);
    Object.values(controls).forEach(el => {
        if (el) el.addEventListener('input', updateUCSCLink)
    });

    // --- Initial Page Setup ---
    if (!hgsid) {
        if(applyLink) applyLink.classList.add('disabled');
        if(resetLink) resetLink.classList.add('disabled');
        const note = document.createElement('p');
        note.innerHTML = '<strong>Note:</strong> Browser view controls are disabled. Open this page from a UCSC track description to enable them.';
        const buttonContainer = document.querySelector('#button-container');
        if (buttonContainer) {
            buttonContainer.insertAdjacentElement('afterend', note);
        }
    } else {
        const resetUrl = new URL(`${ucscHost}/cgi-bin/hgTracks`);
        resetUrl.searchParams.set('hgsid', hgsid);
        resetUrl.searchParams.set(`hg_t_${trackName}_filter`, 'off');
        resetUrl.searchParams.set(trackName, 'dense');
        if(resetLink) resetLink.href = resetUrl.toString();
        updateUCSCLink();
    }
});
