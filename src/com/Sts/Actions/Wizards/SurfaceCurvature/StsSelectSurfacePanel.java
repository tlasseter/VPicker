package com.Sts.Actions.Wizards.SurfaceCurvature;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSelectSurfacePanel extends StsJPanel
{
    private StsSurfaceCurvatureWizard wizard;
    private StsSelectSurface wizardStep;
    public static final String KUWAHARA = "Kuwahara";
    public static String[] FILTERS = {KUWAHARA};
    StsComboBoxFieldBean filterBean;
    String filterString = KUWAHARA;

    private StsModel model = null;

    StsComboBoxFieldBean surfaceBean;
	
    public StsSelectSurfacePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsSurfaceCurvatureWizard)wizard;
    	this.wizardStep = (StsSelectSurface)wizardStep;
    	model  = wizard.getModel();
        buildPanel();
    }

    public void buildPanel()
    {
    	StsSurfaceClass surfaceClass = (StsSurfaceClass)model.getStsClass("com.Sts.DBTypes.StsSurface");
        StsSurface[] surfaces = surfaceClass.getSurfaces();
        StsModelSurfaceClass modelSurfaceClass = (StsModelSurfaceClass)model.getCreateStsClass(StsModelSurface.class);
        StsModelSurface[] modelSurfaces = modelSurfaceClass.getModelSurfaces();

        gbc.fill = gbc.HORIZONTAL;
        gbc.weighty = 0.0;
        if(modelSurfaces.length == 0)
            surfaceBean = new StsComboBoxFieldBean(wizard, "selectedSurface", "Surface:", surfaces);
        else
            surfaceBean = new StsComboBoxFieldBean(wizard, "selectedSurface", "Horizon:", modelSurfaces);

            //surfaceList.initialize(wizard, "selectedSurfaces", "Surfaces:", surfaces);

    	addEndRow(surfaceBean);
		wizard.rebuild();
    }

    public void initialize()
    {
    	surfaceBean.setSelectedIndex(0);
    }
    
}
