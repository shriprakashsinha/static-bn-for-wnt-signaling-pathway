function [uniqueGenes, expressionMatrix, noGenes, noSamples, ...
    groundTruthLabels, transGroundTruthLabels] = readCustomFile(filename)
%
%   Input -
%       filename - name of the file containing the gene expression data
%
%   Output - 
%       uniqueGenes - name of genes gleaned from the file
%       expressionMatrix - 2D matrix containing the gene exp per sample data
%       noGenes - total number of genes available
%       noSamples - total number of samples available
%       groundTruthLabels - original labels available from the files
%       transGroundTruthLabels - labels transformed into numerals
%
%==========================================================================

load(filename)

% extract gene names and remove nan from num matrix
[row, col] = size(custom);
expressionMatrix = [];
uniqueGenes = {};
cnt = 0;
cntDKK = 0;
for i = 1:row
    if ~isnan(num(i,1))
        cnt = cnt + 1;
        if strcmp(custom{i,1}, 'DKK3')
            cntDKK = cntDKK + 1;
            uniqueGenes{cnt} = [custom{i,1}, '-',num2str(cntDKK)];
        else
            uniqueGenes{cnt} = custom{i,1};
        end
        expressionMatrix = [expressionMatrix; num(i,:)];
    end
end

% Find ground truth labels and transform them into numerals
[noGenes, noSamples] = size(expressionMatrix);
groundTruthLabels = {};
transGroundTruthLabels = [];
for i = 1:row
    if ~isempty(strfind(custom{i,1}, 'Sample'))
        j = 1;
        while j <= noSamples
            groundTruthLabels{j} = custom{i,j+1};
            if ~isempty(strfind(groundTruthLabels{j}, 'N'))
                transGroundTruthLabels = [transGroundTruthLabels, -1];
            elseif ~isempty(strfind(groundTruthLabels{j}, 'T'))
                transGroundTruthLabels = [transGroundTruthLabels, +1];
            end
            j = j + 1;
        end
        break
    end
end