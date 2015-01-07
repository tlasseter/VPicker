package com.Sts.UI;

import com.Sts.MVC.StsModel;
import com.Sts.Types.StsSeismicBoundingBox;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.StsParameters;

import javax.swing.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: 3/2/11
 * Time: 9:55 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsSeismicVolume2dGroupBox extends StsGroupBox
{
    StsSeismicBoundingBox volume;
    StsSeismicBoundingBox originalVolume;
    StsFieldBean zDomainBean;
    StsFieldBean sampleTypeBean;
    StsStringFieldBean horizUnitsBean;
    StsStringFieldBean vertUnitsBean;

    StsFloatFieldBean zMinBean;
    StsFloatFieldBean zMaxBean;
    StsFloatFieldBean zIncBean;

    JLabel minLabel = new JLabel("Minimum");
    JLabel maxLabel = new JLabel("Maximum");
    JLabel intervalLabel = new JLabel("Increment");

    JLabel zOrTLabel;

    public StsSeismicVolume2dGroupBox(StsSeismicBoundingBox volume, String title)
    {
  		super(title);
        this.volume = volume;
        buildBeans(false);
        buildPanel();
    }

    public StsSeismicVolume2dGroupBox(StsSeismicBoundingBox volume, String title, StsSeismicBoundingBox originalVolume)
    {
  		super(title);
        this.volume = volume;
        this.originalVolume = originalVolume;
        buildBeans(true);
        buildPanel();
    }

    protected void buildBeans(boolean editable)
    {
        if(!editable)
        {
            zDomainBean = new StsStringFieldBean(volume, "zDomainString", editable, "Domain");
            sampleTypeBean = new StsStringFieldBean(volume, "volumeTypeString", "Sample type");
        }
        else
        {
            String[] zDomainStrings = StsModel.getCurrentProject().getZDomainSupportedStrings();
            zDomainBean = new StsComboBoxFieldBean(this, "zDomainString", zDomainStrings);
            if(StsParameters.volumeTypeIsVelocity(volume.volumeType))
                sampleTypeBean = new StsComboBoxFieldBean(volume, "volumeTypeString", "Sample type", StsParameters.VEL_STRINGS);
            else
                sampleTypeBean = new StsStringFieldBean(volume, "volumeTypeString", "Sample type");
        }
        horizUnitsBean = new StsStringFieldBean(volume, "horzUnitsString", editable, "Horiz Units");
        vertUnitsBean = new StsStringFieldBean(volume, "vertUnitsString", editable, "Vertical Units");

        zMinBean = new StsFloatFieldBean(volume, "zMin", editable);
		zMaxBean = new StsFloatFieldBean(volume, "zMax", editable);
        zIncBean = new StsFloatFieldBean(volume, "zInc", 1, 100*volume.zInc, null, true);
        zIncBean.setEditable(editable);

        if(volume.zDomain.equals(StsParameters.TD_DEPTH_STRING))
            zOrTLabel = new JLabel("Depth:");
        else
            zOrTLabel = new JLabel("Time:");
    }

    protected void buildPanel()
    {
        gbc.fill = HORIZONTAL;
        addUnitsPanel();
        addCoordinatesPanel();
    }

    protected void addUnitsPanel()
    {
        StsJPanel unitsPanel = new StsJPanel();
        unitsPanel.gbc.fill = HORIZONTAL;
        unitsPanel.gbc.anchor = EAST;
        unitsPanel.addToRow(zDomainBean);
        unitsPanel.addEndRow(sampleTypeBean);
        unitsPanel.addToRow(horizUnitsBean);
        unitsPanel.addEndRow(vertUnitsBean);
        add(unitsPanel);
    }

    protected void addCoordinatesPanel()
    {
        StsJPanel coordinatesPanel = new StsJPanel();
        coordinatesPanel.gbc.anchor = CENTER;
		coordinatesPanel.gbc.gridx = 1; // leave empty cell
		coordinatesPanel.addToRow(minLabel);
		coordinatesPanel.addToRow(maxLabel);
		coordinatesPanel.addEndRow(intervalLabel);

        coordinatesPanel.addToRow(zOrTLabel);
		coordinatesPanel.addToRow(zMinBean);
		coordinatesPanel.addToRow(zMaxBean);
		coordinatesPanel.addEndRow(zIncBean);

        add(coordinatesPanel);
    }

    protected String getZDomainString()
    {
        return volume.getZDomainString();
    }

    protected void setZDomainString(String zDomainString)
    {
        String currentZDomainString = volume.getZDomainString();
        if(currentZDomainString == zDomainString) return;
        volume.setZDomainString(zDomainString);
        adjustZRange(zDomainString);
    }

    protected void adjustZRange(String zDomainString)
    {

    }
}
