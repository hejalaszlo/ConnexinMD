function result = calculateDistances(this, respairs)
	result = table.empty;
	chain = unique(this.ResiduePositionCA.chain);
	for iChain = 1:length(chain)
		dist = [];
		for iRespair = 1:length(respairs)
			x1 = this.ResiduePositionCA.x(this.ResiduePositionCA.chain == chain(iChain) & this.ResiduePositionCA.residue == respairs(iRespair,1));
			y1 = this.ResiduePositionCA.y(this.ResiduePositionCA.chain == chain(iChain) & this.ResiduePositionCA.residue == respairs(iRespair,1));
			z1 = this.ResiduePositionCA.z(this.ResiduePositionCA.chain == chain(iChain) & this.ResiduePositionCA.residue == respairs(iRespair,1));

			x2 = this.ResiduePositionCA.x(this.ResiduePositionCA.chain == chain(iChain) & this.ResiduePositionCA.residue == respairs(iRespair,2));
			y2 = this.ResiduePositionCA.y(this.ResiduePositionCA.chain == chain(iChain) & this.ResiduePositionCA.residue == respairs(iRespair,2));
			z2 = this.ResiduePositionCA.z(this.ResiduePositionCA.chain == chain(iChain) & this.ResiduePositionCA.residue == respairs(iRespair,2));

			dist(:,iRespair) = sqrt((x2-x1).^2 + (y2-y1).^2 + (z2-z1).^2);
		end
		result = [result; table((1:size(dist,1))'/10, repmat(chain(iChain), size(dist,1), 1)) array2table(dist)];
	end
	varnames = {'time' 'chain'};
	for iRespair = 1:length(respairs)
		varnames = [varnames strcat('dist_', num2str(respairs(iRespair, 1)), this.Sequence.Res(this.Sequence.Pos == respairs(iRespair, 1)), '-', num2str(respairs(iRespair, 2)), this.Sequence.Res(this.Sequence.Pos == respairs(iRespair, 2)))];
	end
	result.Properties.VariableNames = varnames;
end