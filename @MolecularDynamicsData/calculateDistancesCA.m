function calculateDistancesCA(this)
	this.DistancesCA = table.empty;
	d = table.empty;
	chain = unique(this.ResiduePositionCA.chain);
	res = unique(this.ResiduePositionCA.residue);
	for iRes1 = 1:length(res)
		for iRes2 = iRes1:length(res)
			for iChain1 = 1:length(chain)
				for iChain2 = iChain1:length(chain)
					x1 = this.ResiduePositionCA.x(this.ResiduePositionCA.chain == chain(iChain1) & this.ResiduePositionCA.residue == res(iRes1));
					x2 = this.ResiduePositionCA.x(this.ResiduePositionCA.chain == chain(iChain2) & this.ResiduePositionCA.residue == res(iRes2));
					y1 = this.ResiduePositionCA.y(this.ResiduePositionCA.chain == chain(iChain1) & this.ResiduePositionCA.residue == res(iRes1));
					y2 = this.ResiduePositionCA.y(this.ResiduePositionCA.chain == chain(iChain2) & this.ResiduePositionCA.residue == res(iRes2));
					z1 = this.ResiduePositionCA.z(this.ResiduePositionCA.chain == chain(iChain1) & this.ResiduePositionCA.residue == res(iRes1));
					z2 = this.ResiduePositionCA.z(this.ResiduePositionCA.chain == chain(iChain2) & this.ResiduePositionCA.residue == res(iRes2));
					dist = sqrt((x2-x1).^2 + (y2-y1).^2 + (z2-z1).^2);
					d = [d; table(chain(iChain1), res(iRes1), chain(iChain2), res(iRes2), {dist})];
				 end
			end
		end
	end
	this.DistancesCA = d;
	this.DistancesCA.Properties.VariableNames = {'chain1' 'res1' 'chain2' 'res2' 'dist'};
end 