package com.Sts.UI.Beans;

/**
 * ArrowBeanBeanInfo.java  1.00 97/07/09 Merlin Hughes
 *
 * Copyright (c) 1997 Merlin Hughes, All Rights Reserved.
 *
 * Permission to use, copy, modify, and distribute this software
 * for commercial and non-commercial purposes and without fee is
 * hereby granted provided that this copyright notice appears in
 * all copies.
 *
 * http://prominence.com/                         ego@merlin.org
 */

import java.awt.*;
import java.awt.event.*;
import java.beans.*;

public class ArrowBeanBeanInfo extends SimpleBeanInfo {
  public BeanDescriptor getBeanDescriptor () {
    return new BeanDescriptor (ArrowBean.class, ArrowBeanCustomizer.class);
  }

  public EventSetDescriptor[] getEventSetDescriptors () {
    try {
      EventSetDescriptor[] events = {
        new EventSetDescriptor (ArrowBean.class, "action", ActionListener.class,
                                "actionPerformed")
      };
      return events;
    } catch (IntrospectionException ex) {
      return null;
    }
  }

  public Image getIcon (int iconKind) {
    if (iconKind == ICON_COLOR_16x16) {
      return loadImage ("ArrowBeanIconColor16.gif");
    } else {
      return null;
    }
  }
}
