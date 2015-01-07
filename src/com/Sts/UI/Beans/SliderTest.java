package com.Sts.UI.Beans;

import javax.swing.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Dec 29, 2007
 * Time: 1:16:34 PM
 * To change this template use File | Settings | File Templates.
 */
public class SliderTest
{
    DefaultBoundedRangeModel rangeModel;

    public SliderTest(int min, int max, int value)
    {
        rangeModel = new DefaultBoundedRangeModel(value, 0, min, max);
    }

    public void rangeModelChanged()
    {
        System.out.println("rangeModel changed. Value: " + rangeModel.getValue() + " is adjusting: " + rangeModel.getValueIsAdjusting());
    }

    public Object getRangeModel() { return rangeModel; }
}