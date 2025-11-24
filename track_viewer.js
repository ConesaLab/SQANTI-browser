document.addEventListener('DOMContentLoaded', () => {
    const body = document.body;
    const trackName = body.dataset.trackName;
    let transcripts = [];
    try {
        // Parse the JSON data from the data attribute
        transcripts = JSON.parse(body.dataset.transcripts);
    } catch (e) {
        console.error("Error parsing transcript data:", e);
        // Display an error message to the user
        const container = document.getElementById('container');
        if (container) {
            container.innerHTML = '<h1>Error</h1><p>Could not load transcript data. Please check the data source.</p>';
        }
        return; // Stop execution if data is invalid
    }

    // DOM element references
    const elements = {
        minExons: document.getElementById('min-exons'),
        maxExons: document.getElementById('max-exons'),
        minTss: document.getElementById('min-tss'),
        maxTss: document.getElementById('max-tss'),
        codingStatus: document.getElementById('coding-status'),
        geneSearch: document.getElementById('gene-search'),
        previewButton: document.getElementById('preview-filters'),
        toggleButton: document.getElementById('toggle-table'),
        applyButton: document.getElementById('apply-filters-ucsc'),
        resetButton: document.getElementById('reset-filters-ucsc'),
        resultsCount: document.getElementById('results-count'),
        tableWrapper: document.getElementById('table-wrapper'),
        resultsContainer: document.getElementById('results-container'),
    };

    // --- UTILITY FUNCTIONS ---

    // Get filter values from the input controls
    const getFilterValues = () => ({
        minExons: elements.minExons.value ? parseInt(elements.minExons.value, 10) : null,
        maxExons: elements.maxExons.value ? parseInt(elements.maxExons.value, 10) : null,
        minTss: elements.minTss.value ? parseInt(elements.minTss.value, 10) : null,
        maxTss: elements.maxTss.value ? parseInt(elements.maxTss.value, 10) : null,
        coding: elements.codingStatus.value,
        gene: elements.geneSearch.value.trim().toLowerCase(),
    });

    // Create the HTML table from filtered data
    const createTable = (data) => {
        if (data.length === 0) {
            return '<p>No transcripts match the current filters.</p>';
        }

        const headers = ['Isoform', 'Exons', 'Length', 'Coding', 'Associated Gene', 'TSS Distance', 'TTS Distance'];
        const keys = ['isoform', 'exons', 'length', 'coding', 'associated_gene', 'diff_to_TSS', 'diff_to_TTS'];

        let table = '<table><thead><tr>';
        headers.forEach(header => table += `<th>${header}</th>`);
        table += '</tr></thead><tbody>';

        data.forEach(transcript => {
            table += '<tr>';
            keys.forEach(key => table += `<td>${transcript[key] ?? 'N/A'}</td>`);
            table += '</tr>';
        });

        table += '</tbody></table>';
        return table;
    };

    // --- CORE LOGIC ---

    // Filter transcripts based on current control values
    const filterTranscripts = () => {
        const filters = getFilterValues();
        return transcripts.filter(t => {
            const exonCount = t.exons ?? 0;
            const tssDist = t.diff_to_TSS ?? 0;

            if (filters.minExons !== null && exonCount < filters.minExons) return false;
            if (filters.maxExons !== null && exonCount > filters.maxExons) return false;
            if (filters.minTss !== null && tssDist < filters.minTss) return false;
            if (filters.maxTss !== null && tssDist > filters.maxTss) return false;
            if (filters.coding && filters.coding !== 'all' && t.coding !== filters.coding) return false;
            if (filters.gene && !t.associated_gene?.toLowerCase().includes(filters.gene)) return false;

            return true;
        });
    };

    // Update the results count and table preview
    const updatePreview = () => {
        const filteredData = filterTranscripts();
        elements.resultsCount.textContent = `Showing ${filteredData.length} of ${transcripts.length} transcripts.`;
        elements.resultsContainer.innerHTML = createTable(filteredData);
        updateUcscLink(filteredData.map(t => t.isoform));
    };

    // Update the "Apply to Browser" link with a list of filtered transcript IDs
    const updateUcscLink = (transcriptIds) => {
        if (transcriptIds.length > 0 && transcriptIds.length < transcripts.length) {
            const ucscUrl = `../cgi-bin/hgTracks?db=hg38&hgFind.matches=${transcriptIds.join(',')}`;
            elements.applyButton.href = ucscUrl;
            elements.applyButton.classList.remove('disabled');
        } else {
            elements.applyButton.href = '#';
            elements.applyButton.classList.add('disabled');
        }
    };
    
    // --- EVENT LISTENERS ---

    elements.previewButton.addEventListener('click', updatePreview);

    elements.toggleButton.addEventListener('click', () => {
        const isHidden = elements.tableWrapper.style.display === 'none';
        elements.tableWrapper.style.display = isHidden ? 'block' : 'none';
        if (isHidden) {
            updatePreview(); // Populate table when showing
        }
    });

    // Generate the reset URL
    elements.resetButton.href = `../cgi-bin/hgTracks?db=hg38&setTrackVisibility=hide&${trackName}=dense`;

    // --- INITIALIZATION ---

    // Initial update on page load
    updatePreview(); 
});
