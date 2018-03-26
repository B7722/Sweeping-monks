import QtQuick 2.7
import QtQuick.Controls 2.0

Page {
    width: 600
    height: 400

    title: qsTr("Page 1")

    Button {
        id: button
        x: 390
        y: 302
        text: qsTr("Button")
    }

    TextArea {
        id: textArea
        x: 56
        y: 40
        width: 60
        height: 35
        text: "Items\n"
    }

    TextArea {
        id: textArea1
        x: 56
        y: 103
        text: qsTr("Parameters")
    }

    TextArea {
        id: textArea2
        x: 56
        y: 171
        text: qsTr("Text Area")
    }

    ComboBox {
        id: comboBox
        x: 231
        y: 43
    }

    ComboBox {
        id: comboBox1
        x: 231
        y: 103
    }

    TextField {
        id: textField
        x: 447
        y: 103
        width: 93
        height: 32
        text: qsTr("Text Field")
    }
}
