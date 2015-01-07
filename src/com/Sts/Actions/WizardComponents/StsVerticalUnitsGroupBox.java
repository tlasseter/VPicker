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
public class StsVerticalUnitsGroupBox extends StsGroupBox
{
    public float datumShift = 0.0f;
    public String timeDepthString = StsParameters.TD_DEPTH_STRING;
    public String tvdOrSsString = StsParameters.TVD_STRING;

    public StsVerticalUnitsGroupBox(StsModel model)
    {
        super("Vertical Scale Controls");
        this.timeDepthString = model.getProject().getZDomainString();
        buildPanel();
    }

    private void buildPanel()
    {
        StsComboBoxFieldBean timeDepthBean = new StsComboBoxFieldBean(this, "timeDepthString", "Time or depth:", StsParameters.TD_SELECT_STRINGS);
        timeDepthBean.setToolTipText("Vertical units are either: true depth, seismic depth, or seismic time.  Ignored for binaries.");
        StsComboBoxFieldBean tvdOrSsBean = new StsComboBoxFieldBean(this, "tvdOrSsString", "TVD or SS:", StsParameters.TVD_SUBSEA_STRINGS);
        StsFloatFieldBean datumShiftBean = new StsFloatFieldBean(this, "datumShift", true, "Correction:  ");
        datumShiftBean.setToolTipText("Correction is added to the Depth or Time values");
        datumShiftBean.setColumns(6);

        gbc.fill = HORIZONTAL;
        addBeanToRow(timeDepthBean);
        addBeanToRow(tvdOrSsBean);
        addBeanEndRow(datumShiftBean);
    }

    public float getDatumShift() { return datumShift; }
    public void setDatumShift(float shift) { datumShift = shift; }

    public void setTimeDepthString(String s) { timeDepthString = s; }
    public String getTimeDepthString() { return timeDepthString; }

    public void setTvdOrSsString(String string) { tvdOrSsString = string; }
    public String getTvdOrSsString() { return tvdOrSsString; }
}