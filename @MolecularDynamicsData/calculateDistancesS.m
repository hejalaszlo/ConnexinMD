function calculateDistancesS(this)
	this.DistancesS = table.empty;
	d = table.empty;
	chain = unique(this.ResiduePositionS.chain);
	res = unique(this.ResiduePositionS.residue);
	for iRes1 = 1:length(res)
		for iRes2 = iRes1+1:length(res)
			for iChain = 1:length(chain)
				x1 = this.ResiduePositionS.x(this.ResiduePositionS.chain == chain(iChain) & this.ResiduePositionS.residue == res(iRes1));
				x2 = this.ResiduePositionS.x(this.ResiduePositionS.chain == chain(iChain) & this.ResiduePositionS.residue == res(iRes2));
				y1 = this.ResiduePositionS.y(this.ResiduePositionS.chain == chain(iChain) & this.ResiduePositionS.residue == res(iRes1));
				y2 = this.ResiduePositionS.y(this.ResiduePositionS.chain == chain(iChain) & this.ResiduePositionS.residue == res(iRes2));
				z1 = this.ResiduePositionS.z(this.ResiduePositionS.chain == chain(iChain) & this.ResiduePositionS.residue == res(iRes1));
				z2 = this.ResiduePositionS.z(this.ResiduePositionS.chain == chain(iChain) & this.ResiduePositionS.residue == res(iRes2));
				dist = sqrt((x2-x1).^2 + (y2-y1).^2 + (z2-z1).^2);
				d = [d; table(res(iRes1), res(iRes2), chain(iChain), {dist})];
			end
		end
	end
	this.DistancesS = d;
	this.DistancesS.Properties.VariableNames = {'res1' 'res2' 'chain' 'dist'};
end