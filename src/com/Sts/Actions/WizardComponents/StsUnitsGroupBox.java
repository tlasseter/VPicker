package com.Sts.Actions.WizardComponents;

import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;
import com.Sts.MVC.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jan 22, 2010
 * Time: 9:35:56 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsUnitsGroupBox extends StsGroupBox
{
    public byte vUnits, hUnits, tUnits;

    public StsUnitsGroupBox(StsModel model)
    {
        super("Units");
        hUnits = model.getProject().getXyUnits();
        vUnits = model.getProject().getDepthUnits();
        tUnits = model.getProject().getTimeUnits();
        buildPanel();
    }

    private void buildPanel()
    {
        StsComboBoxFieldBean hUnitsBean = new StsComboBoxFieldBean(this, "horzUnitsString", "Horizontal Units:", StsParameters.DIST_STRINGS);
        StsComboBoxFieldBean vUnitsBean = new StsComboBoxFieldBean(this, "vertUnitsString", "Vertical Units:", StsParameters.DIST_STRINGS);
        StsComboBoxFieldBean tUnitsBean = new StsComboBoxFieldBean(this, "timeUnitsString", "Time Units:", StsParameters.TIME_STRINGS);
        gbc.fill = HORIZONTAL;
        addBeanToRow(hUnitsBean);
        addBeanToRow(vUnitsBean);
        addBeanEndRow(tUnitsBean);
    }

    public String getHorzUnitsString() { return StsParameters.DIST_STRINGS[hUnits]; }
    public String getVertUnitsString() { return StsParameters.DIST_STRINGS[vUnits]; }
    public String getTimeUnitsString() { return StsParameters.TIME_STRINGS[tUnits]; }

    public void setHorzUnitsString(String unitString)
    {
        hUnits = StsParameters.getDistanceUnitsFromString(unitString);
    }

    public void setVertUnitsString(String unitString)
    {
        vUnits = StsParameters.getDistanceUnitsFromString(unitString);
    }

    public void setTimeUnitsString(String unitString)
    {
        tUnits = StsParameters.getTimeUnitsFromString(unitString);
    }

}
