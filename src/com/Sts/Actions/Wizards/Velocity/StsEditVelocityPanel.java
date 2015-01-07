package com.Sts.Actions.Wizards.Velocity;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.SeismicAttributes.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import edu.mines.jtk.interp.SplinesGridder2;
import org.kabeja.dxf.helpers.SplinePoint;

import javax.media.rtp.event.NewParticipantEvent;
import java.awt.*;

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
public class StsEditVelocityPanel extends StsJPanel
{
//    StsFieldBeanPanel editPanel;
    StsGroupBox componentSelectionBox = new StsGroupBox("Select Model Components");
    StsGroupBox datumBox = new StsGroupBox("Define Datums");

    public final static String[] gridStr = {"R-Squared","Spline","Blended","Radial BiHarmonic","v0+kz regression"};
    public final static String[] paramStr={"Radius weight exponent","Spline tension [0 ,, 0.99]","Blended Smoothness [0.25 ... 5]","",""};

    public final static int RSQUARE=0;
    public final static int SPLINES=1;
    public final static int BLENDED=2;
    public final static int RADIAL = 3;
    public final static int v0kz=4;

    public final static int numGridders=5;

    StsListFieldBean surfaceListBean;
    StsComboBoxFieldBean volumeSelectBean;
    StsComboBoxFieldBean gridderComboBoxBean;
    StsBooleanFieldBean useWellControlBean;
    StsBooleanFieldBean fromSonicsBean;
    StsFloatFieldBean topDepthDatumBean;
    StsFloatFieldBean topTimeDatumBean;
	StsFloatFieldBean markerFactorBean;
    StsButton createTestVelocityButton;
    StsModelSurface surface = null;
    StsVelocityWizard wizard = null;
    StsWizardStep wizardStep;
    float topDepthDatum = 0.0f;
    float topTimeDatum = 0.0f;
	float markerFactor = 2.0f;
    int gridType = RSQUARE;
    boolean useWellControl = true;
    boolean useSonics = false;

    double scaleMultiplier = 1.0f;

    static final String noneString = new String("none");

    public StsEditVelocityPanel(StsEditVelocity editVelocity)
    {
		super();
        wizard = (StsVelocityWizard) editVelocity.getWizard();
        wizardStep = editVelocity;
        constructPanel();
    }

    private void constructPanel()
    {
        //fromSonicsBean = new StsBooleanFieldBean(this, "fromSonics", true, "From T/D functions", true);
        //fromSonicsBean.setEditable(false); // Not implemented yet.
        useWellControlBean = new StsBooleanFieldBean(this, "wellControl", true, "Use Well Control:", true);


        gridderComboBoxBean = new StsComboBoxFieldBean(this, "gridTypeString", "Gridder/Extrapolator:", gridStr);
		markerFactorBean = new StsFloatFieldBean(this, "markerFactor", 0.25f, 25.f, "Gridder factor");

        surfaceListBean = new StsListFieldBean(this, "surface", "Horizons:", null);
        createTestVelocityButton = new StsButton("Create Velocity Volume", "Create test velocity volume.", this, "createTestVelocityVolume");
        volumeSelectBean = new StsComboBoxFieldBean(this, "volume", "Velocity volume:", new String[] { noneString });
        topDepthDatumBean = new StsFloatFieldBean(wizard, "topDepthDatum", -10000f, 10000f, "Top depth datum:");
        topTimeDatumBean = new StsFloatFieldBean(wizard, "topTimeDatum", -10000f, 10000f, "Top time datum:");
        componentSelectionBox.gbc.gridwidth = 2;
        componentSelectionBox.gbc.fill = gbc.HORIZONTAL;
        //componentSelectionBox.add(fromSonicsBean);
        componentSelectionBox.add(useWellControlBean);

        componentSelectionBox.add(gridderComboBoxBean);
		componentSelectionBox.add(markerFactorBean);
        componentSelectionBox.gbc.gridwidth = 1;
        componentSelectionBox.add(surfaceListBean);
//        componentSelectionBox.add(createTestVelocityButton);
        componentSelectionBox.add(volumeSelectBean);
        datumBox.add(topTimeDatumBean);
        datumBox.add(topDepthDatumBean);

        this.add(componentSelectionBox,  new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
                ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(3, 2, 0, 0), 0, 0));
        this.add(datumBox,  new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
                ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(3, 2, 0, 0), 0, 0));
    }

    public void initializePanelForSonics()
    {
        useSonics = true;
        setFromSonics(true);
        setWellControl(false);
    }

    public void setSurfaces(StsModelSurface[] surfaces)
    {
        int nSurfaces = 0;
        if (surfaces != null)
        {
            nSurfaces = surfaces.length;
        }

        if (nSurfaces == 0)
        {
            surfaceListBean.setListItems(new Object[] { noneString });
        }
        else
        {
            surfaceListBean.setListItems(surfaces);
            surfaceListBean.setSelectedAll();
        }
    }

    public void setVolumes(StsSeismicVolume[] volumes)
    {
        int nVolumes = 0;
        if (volumes != null)
        {
            nVolumes = volumes.length;
        }
        Object[] volumeObjects = new Object[nVolumes + 1];
        volumeObjects[0] = noneString;
        for (int n = 0; n < nVolumes; n++)
        {
            volumeObjects[n + 1] = volumes[n];
        }
        volumeSelectBean.setListItems(volumeObjects);
        volumeSelectBean.setSelectedIndex(0);
    }

    public StsModelSurface[] getSelectedSurfaces()
    {
        Object[] objects = surfaceListBean.getSelectedObjects();
        if(objects.length == 0)
            return null;
        StsModelSurface[] surfaces = new StsModelSurface[objects.length];
        for(int i=0; i<objects.length; i++)
            surfaces[i] = (StsModelSurface)objects[i];
        return surfaces;
    }

    public void setSurface(Object surface)
    {
        if (surface == noneString)
        {
            this.surface = null;
        }
        else
        {
            this.surface = (StsModelSurface) surface;
        }
    }

	public float getMarkerFactor()
	{
		return this.markerFactor;
	}
	public void setMarkerFactor(float v)
	{
		this.markerFactor = v;
	}
    public Object getSurface()
    {
        return surface;
    }

    public void setFromSonics(boolean value)
    {
        useSonics = value;
        if(useSonics)
        {
            useWellControlBean.setEditable(false);
            surfaceListBean.setEditable(false);
            volumeSelectBean.setEditable(false);
        }
        else
        {
            useWellControlBean.setEditable(true);
            surfaceListBean.setEditable(true);
            volumeSelectBean.setEditable(true);
        }
    }
    public boolean getFromSonics() { return useSonics; }

    public void setWellControl(boolean value)
    {
        useWellControl = value;
        /*
		if (value)
			markerFactorBean.setEditable(true);
		else
			markerFactorBean.setEditable(false);
	    */
    }
    public boolean getWellControl() { return useWellControl; }

    public void setGridType(int value)
        {
            gridType = value;
            markerFactorBean.updateLabel(paramStr[value]);
            markerFactorBean.setText(new String(""));
            markerFactorBean.setVisible(false);
            markerFactorBean.setVisible(true);
            switch (value)   {
                case RSQUARE:
                    markerFactorBean.setText(new String("2")) ;
                    setMarkerFactor(2.f);
                break;
                case SPLINES:
                    markerFactorBean.setText(new String("0.5"));
                    setMarkerFactor(0.5f);
                    break;
                case BLENDED:
                    markerFactorBean.setText(new String(".25"))  ;
                    setMarkerFactor(0.25f);
                    break;
                case v0kz:
                case RADIAL:
                    break;

            }

            /*
            if (value)
                markerFactorBean.setEditable(true);
            else
                markerFactorBean.setEditable(false);
                */
        }

    public int getGridType() { return gridType; }

    public void setGridTypeString(String value)
    {
        for (int i=0; i < numGridders; i++)
        {
            if (value.equals(gridStr[i]))
                setGridType(i);
        }
    }

    public String getGridTypeString()
    {
        return gridStr[getGridType()]  ;
    }

    public void setVolume(Object volume)
    {
        if(volume == noneString)
            wizard.setVelocityVolume(null);
        else
            wizard.setVelocityVolume((StsSeismicVolume) volume);
    }

    public Object getVolume()
    {
        return wizard.getVelocityVolume();
    }

    public void createTestVelocityVolume()
    {
        StsModel model = wizard.getModel();
        StsSeismicVolume seismicVolume = (StsSeismicVolume)model.getStsClass(StsSeismicVolume.class).getFirst();
        StsTestVelocityVolumeConstructor.constructVolume(wizard.getModel(), seismicVolume, true, this);
    }
}
