package com.Sts.Actions.Wizards.Velocity;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;

/**
 * <p>Title: jS2S development</p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author T.J.Lasseter
 * @version c51c
 */
public class StsAssignVelocityPanel extends StsJPanel
{
	StsModel model;
    private StsGroupBox surfaceGroupBox = new StsGroupBox("Surface Definition");
    StsComboBoxFieldBean surfaceListBean;
    StsFloatFieldBean velocityBean;

    private StsGroupBox parmsGroupBox = new StsGroupBox("Model Parameters");
    StsFloatFieldBean zIncrementBean;
    StsFloatFieldBean btmVelocityBean;

	int nSurfaces = 0;
    StsModelSurface[] surfaces = null;
    int selectedSurfaceIndex = -1;
    float selectedSurfaceVelocity = 0.0f;
    float zInc = 4.0f;
    float btmVelocity = 0.0f;
	float[] intervalVelocities; // interval velocities above each surface plus velocity below last surface

    StsVelocityWizard wizard = null;
    StsWizardStep wizardStep;

    static final boolean debug = false;

    public StsAssignVelocityPanel(StsAssignVelocity assignVelocity)
    {
		super();
        wizard = (StsVelocityWizard) assignVelocity.getWizard();
        wizardStep = assignVelocity;
		model = assignVelocity.getModel();
        constructPanel();
    }

    private void constructPanel()
    {
        surfaceListBean = new StsComboBoxFieldBean(this, "surface", "Horizons:");
        surfaceListBean.setToolTipText("Select a horizon");
        velocityBean = new StsFloatFieldBean(this, "velocity", 0f, 20000f, "One-way velocity (" + model.getProject().getVelocityUnits() + ") Above Horizon:");
        velocityBean.setToolTipText("Set velocity above the selected horizon");

        zIncrementBean = new StsFloatFieldBean(this, "zInc", 4.0f, 500.0f, "Model Z Increment:");
        zIncrementBean.setToolTipText("Set the model increment");
        btmVelocityBean = new StsFloatFieldBean(this, "bottomVelocity", 0f, 20000f, "One-way velocity (" + model.getProject().getVelocityUnits() + ") Below Last Horizon:");
        btmVelocityBean.setToolTipText("Set velocity below last horizon");

        surfaceGroupBox.gbc.fill = gbc.HORIZONTAL;
        surfaceGroupBox.addEndRow(surfaceListBean);
        surfaceGroupBox.addEndRow(velocityBean);
        gbc.fill = gbc.HORIZONTAL;
        gbc.gridwidth = 2;
        this.addEndRow(surfaceGroupBox);

        parmsGroupBox.gbc.fill = gbc.HORIZONTAL;
        parmsGroupBox.addEndRow(zIncrementBean);
        parmsGroupBox.addEndRow(btmVelocityBean);
        this.addEndRow(parmsGroupBox);
    }

    public void setSurfaces(StsModelSurface[] surfaces)
    {
        if (surfaces == null || surfaces.length == 0)
		{
			nSurfaces = 0;
			return;
		}
		nSurfaces = surfaces.length;
        this.surfaces = surfaces;
        intervalVelocities = new float[nSurfaces];
        for(int i=0; i<nSurfaces; i++)
            intervalVelocities[i] = 0.0f;
        surfaceListBean.setListItems(surfaces);
    }

    public void setSurface(Object surface)
    {
        if(surfaceListBean != null) selectedSurfaceIndex = surfaceListBean.getSelectedIndex();
        if(debug) System.out.println("setSurface called: selectedSurface set to " + surfaces[selectedSurfaceIndex].getName() + " velocity = " + intervalVelocities[selectedSurfaceIndex]);
        if(velocityBean != null) velocityBean.setValue(intervalVelocities[selectedSurfaceIndex]);
    }

    public Object getSurface()
    {
        if(surfaces == null) return null;
        return surfaces[selectedSurfaceIndex];
    }


    public void setVelocity(float velocity)
    {
		intervalVelocities[selectedSurfaceIndex] = velocity;
        if(debug) System.out.println("setVelocity called: selectedSurface  " + surfaces[selectedSurfaceIndex].getName() + " velocity set to " + velocity);
    }

    public float getVelocity()
    {
        if(intervalVelocities == null) return 0.0f;
        return intervalVelocities[selectedSurfaceIndex];
 //        return selectedSurfaceVelocity;
    }

    public void setZInc(float inc)
    {
        zInc = inc;
    }

    public float getZInc()
    {
        return zInc;
    }
    public void setBottomVelocity(float velocity)
    {
        btmVelocity = velocity;
    }

    public float getBottomVelocity()
    {
        return btmVelocity;
    }

	public float[] getIntervalVelocities()
	{
		return intervalVelocities;
	}
}
