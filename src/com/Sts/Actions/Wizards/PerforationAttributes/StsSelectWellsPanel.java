package com.Sts.Actions.Wizards.PerforationAttributes;

import javax.swing.JScrollPane;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.DBTypes.*;
import com.Sts.MVC.StsModel;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.StsMath;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSelectWellsPanel extends StsJPanel
{
    private StsPerforationAttributesWizard wizard;
    private StsSelectWells wizardStep;

    private StsModel model = null;

    transient StsListFieldBean wellListBean;
	transient JScrollPane wellScrollPane = new JScrollPane();
    transient StsButton clearButton = new StsButton();

    public StsSelectWellsPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsPerforationAttributesWizard)wizard;
    	this.wizardStep = (StsSelectWells)wizardStep;
    	model  = wizard.getModel();
        buildPanel();
    }

    public void buildPanel()
    {
        StsObject[] wells = model.getObjectList(StsWell.class);
        StsWell[] wellsWithPerfs = new StsWell[wells.length];
        int cnt=0;
        for(int i=0; i<wells.length; i++)
        {
            if(((StsWell)wells[i]).getPerforationMarkers() != null)
                wellsWithPerfs[cnt++] = (StsWell)wells[i];
        }
        wellsWithPerfs = (StsWell[])StsMath.trimArray(wellsWithPerfs, cnt);
        wellListBean = new StsListFieldBean(wizard, "selectedWell", null, wellsWithPerfs);
        wellListBean.setSingleSelect();
        wellScrollPane.getViewport().add(wellListBean, null);
        gbc.fill = gbc.BOTH;
        gbc.gridwidth = 2;
        gbc.weighty = 1.0;
		addEndRow(wellScrollPane);
		wizard.rebuild();
    }

    public void initialize()
    {
    	;
    }

}