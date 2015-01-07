
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.FlowSystem.FlowNodes;

import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;

public class StsFlowNodesWizard extends StsWizard
{
	private StsSelectNodeTypes selectTypes;
	private StsRelateNodes relateNodes;

    private StsWizardStep[] mySteps =
    {
        selectTypes = new StsSelectNodeTypes(this),
		relateNodes = new StsRelateNodes(this)
    };

	public StsFlowNodesWizard(StsActionManager actionManager)
    {
        super(actionManager, 500, 400);
        addSteps(mySteps);
    }

    public boolean start()
    {
    	dialog.setTitle("Define Flow Nodes");
    	return super.start();
    }

    public void next()
    {
    	if(currentStep == selectTypes)
			enableFinish();
    	gotoNextStep();
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public boolean end()
    {
		if(!super.end()) return false;
        model.setActionStatus(getClass().getName(), StsModel.STARTED);
		return true;
    }

    static public void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsFlowNodesWizard flowNodesWizard = new StsFlowNodesWizard(actionManager);
        flowNodesWizard.start();
    }
}

