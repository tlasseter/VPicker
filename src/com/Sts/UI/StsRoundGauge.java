
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;


import com.Sts.UI.Beans.StsJPanel;
import com.nextwavesoft.enumeration.*;
import com.nextwavesoft.gauge.*;
import com.nextwavesoft.primitives.GaugeItem;
import com.nextwavesoft.shared.Unit;

import javax.swing.*;
import java.awt.*;

public class StsRoundGauge extends CircularGauge
{
    CircularPointerNeedle speedomainneedle = new CircularPointerNeedle();
    int scale = 1;
    float min, max, highlightMin;
    StsDigitalGauge digitalGauge = null;
    
	public StsRoundGauge()
    {
        this(null, 0, 200, 0, new Color(44,44,44));
    }

	public StsRoundGauge(String title, float min, float max, float highlightMin, Color rimColor)
    {
        super();
        this.min = min;
        this.max = max;
        this.highlightMin = highlightMin;
        configureRoundGauge(title, rimColor);
    }

    private int computeScale()
    {
        // Scale the major tick labels to whole numbers and reflect in label
        if((max - min) < 20)
        {
            while((max - min) < 20)
            {
                max = max * 10;
                min = min * 10;
                if(highlightMin != 0)
                    highlightMin = highlightMin * 10;
                scale *= 10;
            }
            scale = -scale;
        }
        else
        {
            while((max - min) > 2000)
            {
                max = max / 10;
                min = min / 10;
                if(highlightMin != 0)
                    highlightMin = highlightMin / 10;
                scale *= 10;
            }
        }
        return scale;
    }

    public float getMin()
    {
        float minVal = 0.0f;
        if(scale < 0)
            minVal = min / Math.abs(scale);
        else
            minVal = min * Math.abs(scale);
        return minVal; 
    }

    private void configureRoundGauge(String title, Color rimColor)
    {
        scale = computeScale();
        if(scale > 1)
            title = title + " x " + scale;
        else if(scale < 1)
            title = title + " / " + Math.abs(scale);

        //Create circular gauge which holds circular scale.
        setBackground(new Color(44,44,44));
        setRimBrush(rimColor);
        setBackgroundType(CircularBackgroundType.CircularTopGradient);
        setBackgroundRadiusRatio(0.97);

        digitalGauge = new StsDigitalGauge();
        GaugeItem labelItem4 = new GaugeItem(digitalGauge,
		    new Unit(0,UnitType.Percentage),
		    new Unit(20,UnitType.Percentage),
		    new Unit(80,UnitType.Percentage),
		    new Unit(12,UnitType.Percentage));
        addGaugeItem(labelItem4);

        JLabel label1 = new JLabel(title);
        label1.setForeground(Color.white);
        label1.setFont(new Font("Verdana",Font.BOLD,20));
        GaugeItem labelItem = new GaugeItem(label1,
                new Unit(0,UnitType.Percentage),
                new Unit(-20,UnitType.Percentage),
                new Unit(40,UnitType.Percentage),
                new Unit(16,UnitType.Percentage));
        addGaugeItem(labelItem);

        //Create circular scale which holds circular tickset.
        CircularScale circularScale = new CircularScale();
        circularScale.setRadius(new Unit(80,UnitType.Percentage));
        circularScale.setBarExtent(new Unit(2,UnitType.Percentage));
        circularScale.setStartAngle(110);
        circularScale.setSweepAngle(320);
        circularScale.setBarVisible(false);
        getCircularScaleCollection().add(circularScale);

        CircularTickLabelMajor tickLabelMajorMph = new CircularTickLabelMajor();
        tickLabelMajorMph.setFontFamily("Verdana");
        tickLabelMajorMph.setFontSize(new Unit(5,UnitType.Percentage));
        tickLabelMajorMph.setForeground(Color.WHITE);
        tickLabelMajorMph.setOrientation(TextOrientation.Rotated);
        tickLabelMajorMph.setScalePlacement(ScalePlacement.Inside);
        tickLabelMajorMph.setScaleOffset(new Unit(32,UnitType.Percentage));

        //create tick mark major
        CircularTickMarkMajor tickMarkMajorMph = new CircularTickMarkMajor();
        tickMarkMajorMph.setTickMarkAscent(new Unit(3,UnitType.Percentage));
        tickMarkMajorMph.setTickMarkExtent(new Unit(3,UnitType.Percentage));
        tickMarkMajorMph.setBackground(Color.WHITE);
        tickMarkMajorMph.setBorderWidth(0);
        tickMarkMajorMph.setScaleOffset(new Unit(25,UnitType.Percentage));
        tickMarkMajorMph.setScalePlacement(ScalePlacement.Inside);

        CircularTickSet circularTickSet = new CircularTickSet();
        circularTickSet.setMinimum(min);
        circularTickSet.setMaximum(max);

        CircularRange range1 = new CircularRange();
		range1.setStartExtent(new Unit(1,UnitType.Percentage));
		range1.setEndExtent(new Unit(1,UnitType.Percentage));
		range1.setStartValue(min);
		range1.setBackground(Color.green);

        // Colored band around axis
        if(highlightMin != 0)
        {
		    range1.setEndValue(highlightMin);

		    CircularRange range2 = new CircularRange();
		    range2.setStartExtent(new Unit(1,UnitType.Percentage));
		    range2.setEndExtent(new Unit(7,UnitType.Percentage));
		    range2.setStartValue(highlightMin);
		    range2.setEndValue(max);
		    range2.setBackground(Color.red);

		    //add the range to set
		    circularTickSet.getCircularRangeCollection().add(range2);
        }
        else
        {
            range1.setEndValue(max);
        }
        circularTickSet.getCircularRangeCollection().add(range1);

        float minorInterval = (max-min)/100;
        float majorInterval = (max-min)/20;
        circularTickSet.setMinorInterval(minorInterval);
        circularTickSet.setMajorInterval(majorInterval);

        CircularRange circularRange = new CircularRange();
        circularRange.setStartValue(0);
        circularRange.setStartExtent(new Unit(1,UnitType.Percentage));
        circularRange.setEndValue(100);
        circularRange.setEndExtent(new Unit(24,UnitType.Percentage));
        circularRange.setScalePlacement(ScalePlacement.Inside);
        circularRange.setScaleOffset(new Unit(4,UnitType.Percentage));
        circularRange.setBackground(Color.RED);

        //create tick mark minor
        CircularTickMarkMinor tickMarkMinor = new CircularTickMarkMinor();
        tickMarkMinor.setTickMarkAscent(new Unit(0.5,UnitType.Percentage));
        tickMarkMinor.setTickMarkExtent(new Unit(4,UnitType.Percentage));
        tickMarkMinor.setBackground(Color.WHITE);
        tickMarkMinor.setBorderWidth(0);
        tickMarkMinor.setScalePlacement(ScalePlacement.Inside);
        tickMarkMinor.setEndValue(max - 2*minorInterval);

        //create tick mark major
        CircularTickMarkMajor tickMarkMajor = new CircularTickMarkMajor();
        tickMarkMajor.setTickMarkAscent(new Unit(1,UnitType.Percentage));
        tickMarkMajor.setTickMarkExtent(new Unit(6,UnitType.Percentage));
        tickMarkMajor.setBackground(Color.WHITE);
        tickMarkMajor.setBorderWidth(0);
        tickMarkMajor.setScalePlacement(ScalePlacement.Inside);
        //tickMarkMajor.setEndValue(max - 2*minorInterval);

        CircularTickLabelMajor tickLabelMajor = new CircularTickLabelMajor();
        tickLabelMajor.setFontFamily("Verdana");
        tickLabelMajor.setFontSize(new Unit(10,UnitType.Percentage));
        tickLabelMajor.setForeground(Color.WHITE);
        tickLabelMajor.setOrientation(TextOrientation.Normal);
        tickLabelMajor.setScalePlacement(ScalePlacement.Inside);
        tickLabelMajor.setScaleOffset(new Unit(10,UnitType.Percentage));

        //Highlight ticks at extreme
        CircularTickMarkMinor tickMarkMinor2 = null;
        CircularTickMarkMajor tickMarkMajor2 = null;
        if(highlightMin != 0)
        {
            // -----Minor
		    tickMarkMinor2 = new CircularTickMarkMinor();
		    tickMarkMinor2.setTickMarkAscent(new Unit(0.5,UnitType.Percentage));
		    tickMarkMinor2.setTickMarkExtent(new Unit(4,UnitType.Percentage));
		    tickMarkMinor2.setBackground(Color.RED);
		    tickMarkMinor2.setBorderWidth(0);
		    tickMarkMinor2.setScalePlacement(ScalePlacement.Inside);
		    tickMarkMinor2.setStartValue(highlightMin);
		    //tickMarkMinor2.setEndValue(max);

		    // -----Major
		    tickMarkMajor2 = new CircularTickMarkMajor();
		    tickMarkMajor2.setTickMarkAscent(new Unit(1,UnitType.Percentage));
		    tickMarkMajor2.setTickMarkExtent(new Unit(6,UnitType.Percentage));
		    tickMarkMajor2.setBackground(Color.RED);
		    tickMarkMajor2.setBorderWidth(0);
		    tickMarkMajor2.setScalePlacement(ScalePlacement.Inside);
		    tickMarkMajor2.setStartValue(highlightMin);
		    //tickMarkMajor2.setEndValue(max);
        }
        speedomainneedle.setPointerAscent(new Unit(3,UnitType.Percentage));
        speedomainneedle.setPointerExtent(new Unit(105,UnitType.Percentage));

        CircularPointerCap maincap = new CircularPointerCap();
        maincap.setBackground(Color.RED);
        maincap.setPointerExtent(new Unit(25,UnitType.Percentage));
        maincap.setCapType(PointerCapType.CircleConvex);

        circularTickSet.getCircularPointerCollection().add(maincap);
        circularTickSet.getCircularPointerCollection().add(speedomainneedle);

        if(highlightMin != 0)
        {
            circularTickSet.getCircularTickCollection().add(tickMarkMajor2);
            circularTickSet.getCircularTickCollection().add(tickMarkMinor2);
        }
        circularTickSet.getCircularTickCollection().add(tickMarkMajor);
        circularTickSet.getCircularTickCollection().add(tickMarkMinor);
        circularTickSet.getCircularTickCollection().add(tickLabelMajor);

        circularScale.getCircularTickSetCollection().add(circularTickSet);
    }

    public void setValue(float value)
    {
        digitalGauge.setText(Float.toString(value));
        if(scale < 0)
            value = value * Math.abs(scale);
        else
            value = value / scale;
        speedomainneedle.setValue(value);
    }

    public void setDigitalValue(String text)
    {
        digitalGauge.setText(text);
    }

}