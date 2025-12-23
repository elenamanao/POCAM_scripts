#ifndef __PAYLOAD_HEADER_H__
#define __PAYLOAD_HEADER_H__
typedef struct __attribute__((__packed__)) {
  uint16_t nbytes_low;
  uint8_t nbytes_high;
  uint8_t record_type;
  uint32_t icm_ts_low;
  uint16_t icm_ts_high;
  uint8_t encoding_descriptor;
  uint8_t dp_status;
} common_header_t;
#define ICM_EL_ID_STRLEN         2*8U+1
#define PCM_BRD_TEMP_RECORD      0x80
#define PCM_CAL_TRIG_TIMES       0x81
#define PCM_VOLTAGE_DATA         0x82
#define PCM_PD_ADC_DATA          0x83
#define PCM_TDC_DATA             0x84
#define PCM_DEVICEID_DATA        0x85
#define PCM_FPGA_STATE_DATA      0x86
#define PCM_MMB_SENSOR_DATA      0x87
// technical
#define PCM_START_OF_SUBRUN      0x8E
#define PCM_END_OF_FILE          0x8F

typedef struct __attribute__((__packed__)){
  uint16_t subrun;
  uint16_t trigs_per_subrun;
  uint16_t trig_period;
  uint64_t start_time; // ICM clock count
  uint64_t start_epoch_ns;
} pcm_subrun_start_t;

typedef struct __attribute__((__packed__)) {
  uint8_t channel;
  uint8_t status;
  uint16_t raw_value;
  float celcius;
} pcm_brd_temp_record_t;

typedef struct __attribute__((__packed__)) {
  uint8_t source; // encodes 0xTC, where T = target (1, 2) and C = channel (A or B)
  uint16_t n_records;
  uint8_t averaging;
  uint16_t window; // in 100 ns units
} pcm_pd_adc_header_t;

typedef struct __attribute__((__packed__)) {
  uint32_t area;
  uint16_t max_ampl;
  uint16_t base; 
} pcm_pd_adc_record_t;

typedef struct __attribute__((__packed__)){ 
  uint8_t source; // encodes 0xTC, where T = target (1, 2) and C = channel (0...3)
  uint16_t n_records;
  uint8_t mode;
  uint8_t hv;
  uint16_t threshold;
} pcm_tdc_header_t;

typedef struct __attribute__((__packed__)){
  uint8_t channel; // 0..3
  uint8_t ro_err_bits; // 0xRE, where R = epoch roll over, E = decode error
  uint8_t edge; // 0 = rising, 1 = falling
  uint32_t coarse;
  uint8_t fine;
} pcm_tdc_record_t;

typedef struct __attribute__((__packed__)) {
  uint32_t ib_registers[16];
  uint32_t dbm_registers[16];
  uint32_t dbs_registers[16];
  uint32_t ib_pwms[6];
  uint32_t dbm_pwms[4];
  uint32_t dbs_pwms[4];
} pcm_fpga_state_t;

typedef struct __attribute__((__packed__)) {
  char icm_id[ICM_EL_ID_STRLEN];
  char mmb_id[ICM_EL_ID_STRLEN];
  uint16_t icm_fw_version;
  uint64_t pcm_fpga_traceids[3]; // 0=IB, 1=DBM, 2=DBS
  uint16_t pcm_fpga_fw_versions[3]; // 0=IB, 1=DBM, 2=DBS
} pcm_pocam_deviceid_t;

typedef struct __attribute__((__packed__)) {
  uint8_t channel;
  uint16_t raw;
  float voltage;
} pcm_voltage_record_t;

typedef struct __attribute__((__packed__)) {
  float accelerometerXYZ[3];
  float accelerometerT;
  float magnetometerXYZ[3];
  float magnetometerT;
  float pressureP;
  float pressureT;
} pcm_mmb_sensor_record_t;


#endif  // __HAL_MINI_MB_POCAM_DATA_FORMATS_POCAM_RECORD_TYPES_H__
