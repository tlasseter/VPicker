package com.Sts.Actions.Wizards.PreStackExport;

import com.Sts.Actions.Wizards.*;
import com.Sts.UI.*;

public class StsExportHandvel3d extends StsExportHandvel 
{

	public StsExportHandvel3d(StsWizard wizard) 
	{
		super(wizard);
	}

	@Override
	protected StsPreStackHandVelExportPanel createHandVelPanel() 
	{
		return new StsPreStackHandVelExportPanel3d(model, velocityModel);
	}

}
