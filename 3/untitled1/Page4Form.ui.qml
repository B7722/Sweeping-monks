import QtQuick 2.7
import QtQuick.Controls 2.0


Page {
    width: 600
    height: 400

    title: qsTr("Page 4")

    Button {
        id: button
        x: 401
        y: 294
        width: 63
        height: 32
        text: qsTr("Set")
    }

   Label {
        id: textArea
        x: 56
        y: 40
        width: 94
        height: 35
        text: "CAMERA\n"
        verticalAlignment: Text.AlignTop
        horizontalAlignment: Text.AlignHCenter
    }

  Label    {
        id: textArea1
        x: 56
        y: 103
        text: qsTr("Parameters")
    }

   Label{
        id: textArea2
        x: 56
        y: 171
        width: 94
        height: 32
        text: qsTr("Text Area")
    }

    ComboBox {
        id: comboBox
        x: 199
        y: 43
    }

    ComboBox {
        id: comboBox1
        x: 199
        y: 103
    }

    TextField {
        id: textField
        x: 362
        y: 103
        width: 102
        height: 32
    }

    Button {
        id: button1
        x: 486
        y: 294
        width: 63
        height: 32
        text: qsTr("Cancel")
    }
}
