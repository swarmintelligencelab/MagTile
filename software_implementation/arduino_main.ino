/*

The Rutgers Swarm Intelligence Lab

MagTile V1 Master Control Firmware for Arduino Uno

Author: Cy Westbrook (https://westbrook.cy/)
Date: Summer 2024

## Compilation
Designed and tested only on the Arduino UNO R3 platform.
This code should be compiled using the official Arduino IDE (https://www.arduino.cc/en/software)
No additional libraries are needed beyond those included with the Arduino IDE.

## Usage
Connect the SDA and SCL lines of the Arduino UNO to the I2C bus on the MagTile platform.
Connect ~1k pull-up resistors between SDA and 5v, and between SCL and 5v. 
This firmware presents a terminal interface via the Arduino UNO's serial port at 115200 baud.

## Accompanying GUI Configuration Tool
An accompanying web-based GUI Configuration Tool is available:
- Pre-built web page: https://magtile-manager.cy2.me
- Source code: https://github.com/cyfinfaza/sil-magtile-manager-v1
Note: The configuration tool only works in Chromium-based browsers (Chrome, Edge, Brave, Opera, etc.)
as it requires the Web Serial API. (https://caniuse.com/?search=web%20serial%20api)

*/

#include <Arduino.h>
#include <Wire.h>
#include <EEPROM.h>

#define MAX_ADDRESSES 32
#define MAX_ROWS 6
#define MAX_COLUMNS 6
#define POWER_ARRAY_ROWS (MAX_ROWS * 3)
#define POWER_ARRAY_COLUMNS (MAX_COLUMNS * 3)

struct Config {
  uint8_t width;
  uint8_t height;
  uint8_t addressCount;
  uint8_t addressList[MAX_ADDRESSES];
};

Config config;

uint8_t scannedAddresses[MAX_ADDRESSES];
uint8_t scannedAddressCount = 0;

bool blinkingEnabled = false;
unsigned long lastBlinkTime = 0;
bool blinkState = false;

void scanI2CAddresses();
bool testLEDEnable(int address);
bool testLEDDisable(int address);
void setPower(int row, int col, uint16_t power);
void blinkAll();
void enableBlinking();
void disableBlinking();

void setup() {
  Serial.begin(115200);
  Wire.begin();
  pinMode(LED_BUILTIN, OUTPUT);

  // Read config from EEPROM
  EEPROM.get(512, config);

  Serial.println("Command Line Terminal Ready");
}

void loop() {
  if (Serial.available() > 0) {
    String command = Serial.readStringUntil('\n');
    command.trim();  // Remove any trailing newline or whitespace
    processCommand(command);
  }

  if (blinkingEnabled) {
    unsigned long currentTime = millis();
    if (currentTime - lastBlinkTime >= 100) {
      blinkAll();
      lastBlinkTime = currentTime;
    }
  }
}

void processCommand(String command) {
  if (command.startsWith("read_width")) {
    Serial.println("ok : " + String(config.width));
  } else if (command.startsWith("read_height")) {
    Serial.println("ok : " + String(config.height));
  } else if (command.startsWith("write_width")) {
    int value = command.substring(12).toInt();
    if (value >= 0 && value <= MAX_COLUMNS) {
      config.width = value;
      Serial.println("ok : Width set to " + String(config.width));
    } else {
      Serial.println("error");
    }
  } else if (command.startsWith("write_height")) {
    int value = command.substring(13).toInt();
    if (value >= 0 && value <= MAX_ROWS) {
      config.height = value;
      Serial.println("ok : Height set to " + String(config.height));
    } else {
      Serial.println("error");
    }
  } else if (command.startsWith("write_address_list")) {
    String list = command.substring(19);
    config.addressCount = 0;
    int index = 0;
    bool valid = true;
    while (index < list.length()) {
      int spaceIndex = list.indexOf(' ', index);
      if (spaceIndex == -1) {
        spaceIndex = list.length();
      }
      int value = list.substring(index, spaceIndex).toInt();
      if (config.addressCount < MAX_ADDRESSES && value >= 0) {
        config.addressList[config.addressCount++] = value;
      } else {
        valid = false;
        break;
      }
      index = spaceIndex + 1;
    }
    if (valid) {
      Serial.println("ok : Address list updated");
    } else {
      Serial.println("error");
    }
  } else if (command.startsWith("read_address_list")) {
    String response = "ok : ";
    for (int i = 0; i < config.addressCount; i++) {
      response += String(config.addressList[i]);
      if (i < config.addressCount - 1) {
        response += " ";
      }
    }
    Serial.println(response);
  } else if (command.startsWith("scan_addresses")) {
    scanI2CAddresses();
    Serial.print("ok : ");
    for (int i = 0; i < scannedAddressCount; i++) {
      Serial.print(String(scannedAddresses[i], DEC));
      if (i < scannedAddressCount - 1) {
        Serial.print(" ");
      }
    }
    Serial.println();
  } else if (command.startsWith("test_led_enable")) {
    int address = command.substring(16).toInt();
    if (testLEDEnable(address)) {
      Serial.println("ok : LED enabled at address " + String(address));
    } else {
      Serial.println("error");
    }
  } else if (command.startsWith("test_led_disable")) {
    int address = command.substring(17).toInt();
    if (testLEDDisable(address)) {
      Serial.println("ok : LED disabled at address " + String(address));
    } else {
      Serial.println("error");
    }
  } else if (command.startsWith("store_config")) {
    EEPROM.put(512, config);
    Serial.println("ok : Configuration stored to EEPROM");
  } else if (command.startsWith("set_power")) {
    int firstSpace = command.indexOf(' ', 10);
    int secondSpace = command.indexOf(' ', firstSpace + 1);
    int row = command.substring(10, firstSpace).toInt();
    int col = command.substring(firstSpace + 1, secondSpace).toInt();
    int power = command.substring(secondSpace + 1).toInt();
    if (row >= 0 && row < POWER_ARRAY_ROWS && col >= 0 && col < POWER_ARRAY_COLUMNS && power >= 0 && power <= 4096) {
      setPower(row, col, power);
      Serial.println("ok : Power set at (" + String(row) + "," + String(col) + ") to " + String(power));
    } else {
      Serial.println("error");
    }
  } else if (command.startsWith("get_power")) {
    Serial.println("error");
  } else if (command == "blinkall_start") {
    enableBlinking();
    Serial.println("ok : Blinking started");
  } else if (command == "blinkall_stop") {
    disableBlinking();
    Serial.println("ok : Blinking stopped");
  } else {
    Serial.println("error");
  }
}

void scanI2CAddresses() {
  scannedAddressCount = 0;
  for (int address = 1; address <= 127; address++) {
    Wire.beginTransmission(address);
    if (Wire.endTransmission() == 0 && address != 112) {
      scannedAddresses[scannedAddressCount++] = address;

      Wire.beginTransmission(address);
      Wire.write(0x00);
      Wire.write(B00110001);  // MODE 1
      Wire.endTransmission();

      Wire.beginTransmission(address); 
      Wire.write(0xFE);
      Wire.write(8);  // Timer Prescaler (3-255)
      Wire.endTransmission();

      Wire.beginTransmission(address);
      Wire.write(0x00);
      Wire.write(B10100001);  // MODE 1
      Wire.endTransmission();

      Wire.beginTransmission(address);
      Wire.write(0x01);
      Wire.write(B00000100);  // MODE 2
      Wire.endTransmission();
    }
  }
}

bool testLEDEnable(int address) {
  Wire.beginTransmission(address);
  Wire.write(0x2A);
  Wire.write(0x00);
  Wire.write(0x00);
  Wire.write(0xFF); // On LS2B
  Wire.write(0x0F); // On MSB
  Wire.endTransmission();
  return true;
}

bool testLEDDisable(int address) {
  Wire.beginTransmission(address);
  Wire.write(0x2A);
  Wire.write(0x00);
  Wire.write(0x00);
  Wire.write(0x00); // On LS2B
  Wire.write(0x00); // On MSB
  Wire.endTransmission();
  return true;
}

void setPower(int row, int col, uint16_t power) {
  // powerLevels[row][col] = power;
  uint8_t tileRow = row/3;
  uint8_t tileCol = col/3;
  uint8_t rowOnTile = row - tileRow*3;
  uint8_t colOnTile = col - tileCol*3;
  uint8_t tileAddr = config.addressList[tileRow*config.width+tileCol];
  uint8_t coilAddr = 0x06 + 4 * (rowOnTile * 3 + colOnTile);
  // Serial.println(tileRow);
  // Serial.println(tileCol);
  // Serial.println(rowOnTile);
  // Serial.println(colOnTile);
  // Serial.println(tileAddr);
  // Serial.println(coilAddr);
  // Serial.println(power & 0x00FF);
  // Serial.println((power & 0x0F00) >> 8);
  Wire.beginTransmission(tileAddr);
  Wire.write(coilAddr);
  Wire.write(0x00);
  Wire.write(0x00);
  Wire.write(power & 0x00FF); // On LS2B
  Wire.write((power & 0x0F00) >> 8); // On MSB
  Wire.endTransmission();
}

void blinkAll() {
  blinkState = !blinkState;
  // Serial.println("Blinking: " + String(blinkState ? "ON" : "OFF"));
  for (int i = 0; i < scannedAddressCount; i++) {
    Wire.beginTransmission(scannedAddresses[i]);
    Wire.write(0xFA);  // ALL_LED_ON_L register
    Wire.write(0x00);
    Wire.write(0x00);
    if (blinkState) {
      Wire.write(0xFF);
      Wire.write(0x0F);  // Set all LEDs fully on
    } else {
      Wire.write(0x00);
      Wire.write(0x00);  // Set all LEDs fully off
    }
    Wire.endTransmission();
  }
}

void enableBlinking() {
  blinkingEnabled = true;
  lastBlinkTime = millis();
}

void disableBlinking() {
  blinkingEnabled = false;
  // Turn off all LEDs when stopping the blink
  for (int i = 0; i < scannedAddressCount; i++) {
    Wire.beginTransmission(scannedAddresses[i]);
    Wire.write(0xFA);  // ALL_LED_ON_L register
    Wire.write(0x00);
    Wire.write(0x00);
    Wire.write(0x00);
    Wire.write(0x00);  // Set all LEDs fully off
    Wire.endTransmission();
  }
}
