function pres_corr = PeGS_SortByID(pres_corr)       %% sorting pres_corr-structure on ID (column 1)
        pres_corr_fields = fieldnames(pres_corr);
        pres_corr_cell = struct2cell(pres_corr);
        sz = size(pres_corr_cell);            % Notice that the this is a 3 dimensional array.
        pres_corr_cell = reshape(pres_corr_cell, sz(1), []);      % Px(MxN)
        pres_corr_cell = pres_corr_cell';                         % (MxN)xP
        pres_corr_cell = sortrows(pres_corr_cell, 1);
        pres_corr_cell = reshape(pres_corr_cell', sz);
        pres_corr = cell2struct(pres_corr_cell, pres_corr_fields, 1);
end