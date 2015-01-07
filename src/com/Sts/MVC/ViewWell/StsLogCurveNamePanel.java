package com.Sts.MVC.ViewWell;

import com.Sts.DBTypes.StsColor;
import com.Sts.DBTypes.StsLogCurve;
import com.Sts.DBTypes.StsLogCurveType;
import com.Sts.UI.Beans.StsJPanel;
import com.Sts.Utilities.StsMath;

import javax.swing.*;
import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: 2/22/11
 * Time: 11:56 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsLogCurveNamePanel extends StsJPanel
{
    public StsLogCurveType logCurveType;
    public JRadioButton selectCurveRadioButton = new JRadioButton();
    private JLabel minValueLabel = new JLabel();
    private JLabel curveLabel = new JLabel();
    private JLabel maxValueLabel = new JLabel();
    private JLabel scaleTypeLabel = new JLabel();

    public StsLogCurveNamePanel(StsLogCurve logCurve)
    {
        super();
        logCurveType = logCurve.getLogCurveType();
        String curveName = logCurve.getName();
        setName(curveName);
        setBorder(BorderFactory.createRaisedBevelBorder());
        //        setMinimumSize(new Dimension(curveTrackWidth, 25));
        //        setPreferredSize(new Dimension(curveTrackWidth, 25));
        selectCurveRadioButton.setMargin(new Insets(0, 0, 0, 2));
//        selectCurveRadioButton.setPreferredSize(new Dimension(15, 15));
//        selectCurveRadioButton.setMinimumSize(new Dimension(15, 15));
        selectCurveRadioButton.setContentAreaFilled(false);
        // track.addRadioButtonToGroup(selectCurveRadioButton);
        selectCurveRadioButton.setActionCommand(curveName);

        curveLabel.setText(curveName);
        if(logCurveType != null)
        {
            setValueLabels();
            StsColor color = logCurve.getStsColor();
            setBackground(color.getColor());
            selectCurveRadioButton.setForeground(color.getColor());
        }

        addToRow(selectCurveRadioButton);
        gbc.fill = GridBagConstraints.HORIZONTAL;
        addToRow(minValueLabel);
        gbc.fill = GridBagConstraints.NONE;
        addToRow(curveLabel);
        gbc.fill = GridBagConstraints.HORIZONTAL;
        addToRow(maxValueLabel);
        gbc.fill = GridBagConstraints.NONE;
        addToRow(scaleTypeLabel);
    }

    public StsLogCurveNamePanel(String curveName, float min, float max, String scaleType)
    {
        super();
        setName(curveName);
        setBorder(BorderFactory.createRaisedBevelBorder());
        //        setMinimumSize(new Dimension(curveTrackWidth, 25));
        //        setPreferredSize(new Dimension(curveTrackWidth, 25));
        selectCurveRadioButton.setMargin(new Insets(0, 0, 0, 2));
//        selectCurveRadioButton.setPreferredSize(new Dimension(15, 15));
//        selectCurveRadioButton.setMinimumSize(new Dimension(15, 15));
        selectCurveRadioButton.setContentAreaFilled(false);
        // track.addRadioButtonToGroup(selectCurveRadioButton);
        selectCurveRadioButton.setActionCommand(curveName);

        curveLabel.setText(curveName);
        setValueLabels();
        scaleTypeLabel.setText(scaleType);
        Color color = Color.WHITE;
        setBackground(color);

        addToRow(selectCurveRadioButton);
        gbc.fill = GridBagConstraints.HORIZONTAL;
        addToRow(minValueLabel);
        gbc.fill = GridBagConstraints.NONE;
        addToRow(curveLabel);
        gbc.fill = GridBagConstraints.HORIZONTAL;
        addToRow(maxValueLabel);
        gbc.fill = GridBagConstraints.NONE;
        addToRow(scaleTypeLabel);
    }

    public void setValueLabels()
    {
        String minValueString = StsMath.formatNumber(logCurveType.getScaleMin(), 5, 5);
        minValueLabel.setText(minValueString);
        String maxValueString = StsMath.formatNumber(logCurveType.getScaleMax(), 5, 5);
        maxValueLabel.setText(maxValueString);
        scaleTypeLabel.setText(logCurveType.getScaleTypeString());
    }

    public void setMinValueLabel(String minValueLabel)
    {
        this.minValueLabel.setText(minValueLabel);
    }

    public void setMaxValueLabel(String maxValueLabel)
    {
        this.maxValueLabel.setText(maxValueLabel);
    }

    public void setScaleTypeLabel(String scaleTypeLabel)
    {
        this.scaleTypeLabel.setText(scaleTypeLabel);
    }
}
