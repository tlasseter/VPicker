package com.Sts.UI.Beans;

import java.beans.*;

/**
 * <p>Title: Field Beans Development</p>
 * <p>Description: General beans for generic panels.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.0
 */

public class StsRadioButtonFieldBeanBeanInfo extends SimpleBeanInfo
{
    private Class beanClass = StsRadioButtonFieldBean.class;
    private String iconColor16x16Filename;
    private String iconColor32x32Filename = "StsRadioButtonFieldBean.gif";
    private String iconMono16x16Filename;
    private String iconMono32x32Filename;

    public StsRadioButtonFieldBeanBeanInfo()
    {
    }
    public PropertyDescriptor[] getPropertyDescriptors()
    {
        try
        {
            PropertyDescriptor _horizontalAlignment = new PropertyDescriptor("horizontalAlignment", beanClass, null, "setHorizontalAlignment");
            PropertyDescriptor _horizontalTextPosition = new PropertyDescriptor("horizontalTextPosition", beanClass, null, "setHorizontalTextPosition");
            PropertyDescriptor _selected = new PropertyDescriptor("selected", beanClass, "isSelected", "setSelected");
            PropertyDescriptor _text = new PropertyDescriptor("text", beanClass, "getText", "setText");
            PropertyDescriptor[] pds = new PropertyDescriptor[] {
	            _horizontalAlignment,
	            _horizontalTextPosition,
	            _selected,
	            _text};
            return pds;
        }
        catch(IntrospectionException ex)
        {
            ex.printStackTrace();
            return null;
        }
    }
    public java.awt.Image getIcon(int iconKind)
    {
        switch (iconKind) {
            case BeanInfo.ICON_COLOR_16x16:
              return iconColor16x16Filename != null ? loadImage(iconColor16x16Filename) : null;
            case BeanInfo.ICON_COLOR_32x32:
              return iconColor32x32Filename != null ? loadImage(iconColor32x32Filename) : null;
            case BeanInfo.ICON_MONO_16x16:
              return iconMono16x16Filename != null ? loadImage(iconMono16x16Filename) : null;
            case BeanInfo.ICON_MONO_32x32:
              return iconMono32x32Filename != null ? loadImage(iconMono32x32Filename) : null;
        }
        return null;
    }
    public BeanInfo[] getAdditionalBeanInfo()
    {
        Class superclass = beanClass.getSuperclass();
        try
        {
            BeanInfo superBeanInfo = Introspector.getBeanInfo(superclass);
            return new BeanInfo[] { superBeanInfo };
        }
        catch(IntrospectionException ex)
        {
            ex.printStackTrace();
            return null;
        }
    }
}